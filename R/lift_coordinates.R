#' Lift m6A Transcript Coordinates to Genomic Coordinates
#'
#' Converts transcript-level m6A coordinates (as reported by m6Anet)
#' to genomic coordinates using a GTF/GFF reference. This is the critical step
#' required to compare m6A sites across isoforms: transcript position 245 in
#' ENST001 and position 245 in ENST002 refer to different genomic loci, so
#' comparisons must be made in genomic space.
#'
#' @param m6a_sites data.table with columns: transcript_id, position,
#'   probability (and any additional columns from parse_m6anet())
#' @param gtf_file Path to GTF or GFF file for the reference genome/transcriptome
#'
#' @return GRanges object with one range per successfully mapped m6A site.
#'   Metadata columns include:
#'   \describe{
#'     \item{transcript_id}{Original transcript identifier}
#'     \item{transcript_position}{Original 1-indexed position in transcript}
#'     \item{probability}{m6A modification probability}
#'   }
#'   Any additional columns from \code{m6a_sites} are forwarded as metadata.
#'
#' @details
#' The GTF file specifies which exons belong to each transcript and their
#' positions on the chromosome. This function uses that map to convert from
#' transcript numbering to chromosome numbering.
#'
#' mapFromTranscripts() requires the exon structure to correctly skip introns.
#' Using transcript boundaries alone would place position 150 inside an intron:
#'
#' \preformatted{
#'   Isoform A exon structure (from GTF):
#'     Exon1: chr10:1000-1100  (transcript positions 1-100)
#'     Exon2: chr10:5000-5200  (transcript positions 101-300)
#'     [Intron chr10:1101-4999 - not part of the transcript]
#'
#'   Transcript position 150 = 50 bases into Exon2
#'     = chr10:5000 + 50 - 1 = chr10:5049   CORRECT
#'
#'   Using transcript boundaries only:
#'     = 1000 + 150 = chr10:1150            WRONG (inside intron)
#' }
#'
#' @section Enabling condition-level comparison:
#' To compare conditions downstream, add a \code{condition} column to each
#' parsed m6A table and combine them \emph{without} deduplicating:
#' \preformatted{
#' m6a_a[, condition := "WT"]
#' m6a_b[, condition := "MUT"]
#' combined <- rbind(m6a_a, m6a_b)
#' gr <- lift_m6a_to_genomic(combined, gtf)
#' }
#' Extra columns are forwarded automatically. Without a \code{condition}
#' column, \code{annotate_m6a_switches_genomic()} can only report isoform
#' membership, not condition-specific methylation.
#'
#' @examples
#' \dontrun{
#' m6a_sites <- parse_m6anet("m6anet_predictions.csv")
#' genomic_sites <- lift_m6a_to_genomic(m6a_sites, "genome.gtf")
#' }
#'
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicFeatures exonsBy mapFromTranscripts
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
#' @importFrom S4Vectors mcols
#' @import data.table
#'
#' @export
lift_m6a_to_genomic <- function(m6a_sites, gtf_file) {

  if (!data.table::is.data.table(m6a_sites)) {
    stop("m6a_sites must be a data.table from parse_m6anet()")
  }
  required_cols <- c("transcript_id", "position", "probability")
  if (!all(required_cols %in% names(m6a_sites))) {
    stop("m6a_sites must contain columns: ", paste(required_cols, collapse = ", "))
  }
  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }

  # Build TxDb from GTF
  txdb <- txdbmaker::makeTxDbFromGFF(gtf_file)

  # Exon structure per transcript. mapFromTranscripts() needs this to know
  # where introns are; transcript boundaries alone are not sufficient.
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

  other_cols <- setdiff(names(m6a_sites), c("transcript_id", "position", "probability"))

  result_list <- list()
  missing_tx  <- data.table::data.table(transcript_id = character(), n = integer())
  failed_map  <- data.table::data.table(transcript_id = character(), n = integer())
  error_msgs  <- character()

  sites_by_tx <- split(m6a_sites, by = "transcript_id", keep.by = FALSE, drop = TRUE)

  for (tx_id in names(sites_by_tx)) {
    tx_sites <- sites_by_tx[[tx_id]]

    # NOTE: single bracket. `all_exons[[tx_id]]` returns a plain GRanges with
    # names dropped, and mapFromTranscripts() errors with
    # "'transcripts' must have names". Single bracket keeps a length-1
    # GRangesList with the name intact.
    tx_range <- all_exons[tx_id]

    if (is.null(tx_range) || length(tx_range) == 0) {
      missing_tx <- rbind(
        missing_tx,
        data.table::data.table(transcript_id = tx_id, n = nrow(tx_sites))
      )
      next
    }

    tx_pos    <- as.integer(tx_sites$position)
    valid_pos <- !is.na(tx_pos)
    if (!all(valid_pos)) {
      failed_map <- rbind(
        failed_map,
        data.table::data.table(transcript_id = tx_id, n = sum(!valid_pos))
      )
    }
    if (!any(valid_pos)) next

    tx_sites_valid <- tx_sites[valid_pos]
    tx_pos_valid   <- tx_pos[valid_pos]

    tx_query <- GenomicRanges::GRanges(
      seqnames = rep(tx_id, length(tx_pos_valid)),
      ranges   = IRanges::IRanges(start = tx_pos_valid, end = tx_pos_valid)
    )

    genomic_ranges <- tryCatch(
      GenomicFeatures::mapFromTranscripts(tx_query, tx_range),
      error = function(e) {
        error_msgs <<- c(error_msgs, conditionMessage(e))
        failed_map <<- rbind(
          failed_map,
          data.table::data.table(transcript_id = tx_id, n = nrow(tx_sites_valid))
        )
        NULL
      }
    )

    if (is.null(genomic_ranges) || length(genomic_ranges) == 0) {
      failed_map <- rbind(
        failed_map,
        data.table::data.table(transcript_id = tx_id, n = nrow(tx_sites_valid))
      )
      next
    }

    query_hits <- if ("xHits" %in% names(S4Vectors::mcols(genomic_ranges))) {
      as.integer(S4Vectors::mcols(genomic_ranges)$xHits)
    } else {
      seq_len(min(length(genomic_ranges), nrow(tx_sites_valid)))
    }

    mapped_sites <- tx_sites_valid[query_hits]
    if (nrow(mapped_sites) < nrow(tx_sites_valid)) {
      failed_map <- rbind(
        failed_map,
        data.table::data.table(transcript_id = tx_id,
                               n = nrow(tx_sites_valid) - nrow(mapped_sites))
      )
    }

    genomic_ranges$transcript_id       <- mapped_sites$transcript_id
    genomic_ranges$transcript_position <- as.integer(mapped_sites$position)
    genomic_ranges$probability         <- mapped_sites$probability

    for (col in other_cols) {
      genomic_ranges[[col]] <- mapped_sites[[col]]
    }

    result_list[[length(result_list) + 1L]] <- genomic_ranges
  }

  result_list <- Filter(Negate(is.null), result_list)

  if (nrow(missing_tx) > 0) {
    miss_sum <- missing_tx[, .(n = sum(n)), by = transcript_id]
    warning(sprintf(
      "Skipped %d site(s): %d transcript(s) not found in GTF (%s).",
      sum(miss_sum$n), nrow(miss_sum),
      paste(utils::head(miss_sum$transcript_id, 5), collapse = ", ")
    ), call. = FALSE)
  }

  if (nrow(failed_map) > 0) {
    fail_sum <- failed_map[, .(n = sum(n)), by = transcript_id]
    warning(sprintf(
      "Could not map %d site(s) across %d transcript(s).",
      sum(fail_sum$n), nrow(fail_sum)
    ), call. = FALSE)
  }

  # Surface the actual error messages, not just a count. A silent tryCatch
  # previously hid a systematic failure behind an aggregate warning.
  if (length(error_msgs) > 0) {
    warning(sprintf(
      "mapFromTranscripts() raised %d error(s). Distinct messages:\n  %s",
      length(error_msgs),
      paste(unique(error_msgs), collapse = "\n  ")
    ), call. = FALSE)
  }

  if (length(result_list) == 0) {
    stop(paste(
      "No m6A sites could be mapped to genomic coordinates.",
      "Check GTF file and transcript IDs."
    ))
  }

  do.call(c, result_list)
}


#' Annotate Isoform Switches Using Genomic m6A Coordinates
#'
#' Compares m6A sites at the same genomic locus across the two isoforms of a
#' switch pair, and - when condition information is supplied - between the two
#' experimental conditions. These are two distinct questions and are reported
#' in two separate columns.
#'
#' @param m6a_sites_gr GRanges object from \code{lift_m6a_to_genomic()}.
#'   Must have metadata columns \code{transcript_id} and \code{probability}.
#'   If a \code{condition} column is present, condition-level fates are also
#'   computed.
#' @param iso_switches data.table from \code{parse_isoform_switch()}.
#'   Should include condition_1, condition_2, and dif columns.
#'
#' @return data.table with two independent classifications:
#'   \describe{
#'     \item{isoform_status}{STRUCTURAL. Which transcript of the pair carries
#'       the site: \code{ISOFORM_A_ONLY}, \code{ISOFORM_B_ONLY},
#'       \code{IN_BOTH_ISOFORMS}. Always computed.}
#'     \item{m6a_fate}{REGULATORY. Which condition detected the site:
#'       \code{LOST} (condition_1 only), \code{GAINED} (condition_2 only),
#'       \code{RETAINED} (both). NA unless a \code{condition} column was
#'       supplied to \code{lift_m6a_to_genomic()}.}
#'     \item{probability_a, probability_b}{m6A probability on each isoform}
#'     \item{probability_c1, probability_c2}{m6A probability in each condition}
#'   }
#'   Plus: gene_id, isoform_a, isoform_b, condition_1, condition_2,
#'   genomic_position, seqname, start, end, transcript_position, strand,
#'   m6a_in_isoform_a, m6a_in_isoform_b, fdr, dif.
#'
#' @section Two different questions:
#' \code{isoform_status} answers "do the two isoforms of this switch carry
#' different m6A sites?" This is a consequence of transcript structure - a
#' longer or differently spliced isoform will carry more sites.
#'
#' \code{m6a_fate} answers "does methylation at this position differ between
#' conditions?" This is a regulatory question, independent of isoform usage.
#'
#' These are easily conflated. A site can be detected in both conditions while
#' sitting on only one isoform of the pair - in which case
#' \code{isoform_status} is ISOFORM_B_ONLY but \code{m6a_fate} is RETAINED.
#'
#' @section Supplying condition information:
#' \preformatted{
#' m6a_a[, condition := "WT"]
#' m6a_b[, condition := "MUT"]
#' combined <- rbind(m6a_a, m6a_b)   # do NOT deduplicate across conditions
#' gr  <- lift_m6a_to_genomic(combined, gtf)
#' res <- annotate_m6a_switches_genomic(gr, iso_switches)
#' }
#' Condition labels must match \code{condition_1} / \code{condition_2} in
#' \code{iso_switches}. A warning is issued if they do not.
#'
#' @section Why genomic coordinates:
#' Transcript position 150 in two isoforms with different exon structures maps
#' to different genomic locations. Comparing transcript positions directly
#' produces false calls:
#' \preformatted{
#'   Isoform A: Exon1 (chr10:1000-1100) + Exon2 (chr10:5000-5200)
#'   Isoform B: Exon1 (chr10:1000-1100) + Exon3 (chr10:8000-8200)
#'
#'   By transcript position:  150 in A has m6A, 150 in B has m6A
#'                            -> "both"            WRONG
#'
#'   By genomic coordinate:   chr10:5049 in A has m6A
#'                            chr10:5049 absent from B
#'                            -> ISOFORM_A_ONLY    CORRECT
#' }
#'
#' @examples
#' \dontrun{
#' m6a_a <- parse_m6anet("cond_a.csv"); m6a_a[, condition := "WT"]
#' m6a_b <- parse_m6anet("cond_b.csv"); m6a_b[, condition := "MUT"]
#' iso   <- parse_isoform_switch("switches.txt")
#'
#' gr  <- lift_m6a_to_genomic(rbind(m6a_a, m6a_b), "genome.gtf")
#' res <- annotate_m6a_switches_genomic(gr, iso)
#'
#' table(res$isoform_status)                    # structural
#' table(res$m6a_fate, useNA = "ifany")         # regulatory
#' table(res$isoform_status, res$m6a_fate)      # are they independent?
#' }
#'
#' @importFrom GenomicRanges findOverlaps seqnames start end strand GRanges
#' @importFrom S4Vectors subjectHits mcols
#' @import data.table
#' @import methods
#'
#' @export
annotate_m6a_switches_genomic <- function(m6a_sites_gr, iso_switches) {

  if (!methods::is(m6a_sites_gr, "GRanges")) {
    stop("m6a_sites_gr must be a GRanges object from lift_m6a_to_genomic()")
  }
  if (!"transcript_id" %in% names(S4Vectors::mcols(m6a_sites_gr))) {
    stop("m6a_sites_gr must have a 'transcript_id' metadata column")
  }
  if (!"probability" %in% names(S4Vectors::mcols(m6a_sites_gr))) {
    stop("m6a_sites_gr must have a 'probability' metadata column")
  }
  if (!data.table::is.data.table(iso_switches)) {
    stop("iso_switches must be a data.table from parse_isoform_switch()")
  }

  if (nrow(iso_switches) == 0) {
    warning("iso_switches is empty. Returning empty result.", call. = FALSE)
    return(data.table::data.table())
  }

  # ── Is condition information available? ─────────────────────────────────────
  has_condition <- "condition" %in% names(S4Vectors::mcols(m6a_sites_gr))

  if (!has_condition) {
    warning(
      "No 'condition' column found in m6a_sites_gr.\n",
      "  isoform_status will be computed (structural: which isoform carries the site).\n",
      "  m6a_fate will be NA (regulatory: which condition detected the site).\n",
      "  To enable condition comparison, add a 'condition' column to each parsed\n",
      "  m6A table, rbind them WITHOUT deduplicating, and lift the combined table.",
      call. = FALSE
    )
  } else {
    obs_conds <- unique(m6a_sites_gr$condition)
    exp_conds <- unique(c(iso_switches$condition_1, iso_switches$condition_2))
    exp_conds <- exp_conds[!is.na(exp_conds)]
    if (length(exp_conds) > 0 && !all(exp_conds %in% obs_conds)) {
      warning("Condition labels in iso_switches (",
              paste(exp_conds, collapse = ", "),
              ") do not all appear in m6a_sites_gr (",
              paste(obs_conds, collapse = ", "),
              "). m6a_fate may be incorrect.", call. = FALSE)
    }
  }

  m6a_by_tx   <- split(m6a_sites_gr, m6a_sites_gr$transcript_id)
  result_list <- list()

  for (i in seq_len(nrow(iso_switches))) {
    switch_row <- iso_switches[i]
    gene  <- switch_row$gene_id
    iso_a <- switch_row$isoform_a
    iso_b <- switch_row$isoform_b
    fdr   <- switch_row$fdr

    condition_1 <- if ("condition_1" %in% names(switch_row))
      switch_row$condition_1 else NA_character_
    condition_2 <- if ("condition_2" %in% names(switch_row))
      switch_row$condition_2 else NA_character_
    dif <- if ("dif" %in% names(switch_row)) switch_row$dif else NA_real_

    m6a_in_a_gr <- m6a_by_tx[[iso_a]]
    m6a_in_b_gr <- m6a_by_tx[[iso_b]]
    if (is.null(m6a_in_a_gr)) m6a_in_a_gr <- GenomicRanges::GRanges()
    if (is.null(m6a_in_b_gr)) m6a_in_b_gr <- GenomicRanges::GRanges()

    if (length(m6a_in_a_gr) == 0 && length(m6a_in_b_gr) == 0) next

    # Distinct single-base genomic loci; do not merge adjacent sites
    merged <- unique(c(m6a_in_a_gr, m6a_in_b_gr))

    for (j in seq_len(length(merged))) {
      genomic_pos <- merged[j]

      # ── STRUCTURAL: which isoform of the pair carries the site ──────────────
      overlaps_a <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_a_gr, type = "equal")
      in_a       <- length(overlaps_a) > 0
      prob_a     <- if (in_a) {
        m6a_in_a_gr$probability[S4Vectors::subjectHits(overlaps_a)[1]]
      } else NA_real_
      tx_pos_a   <- if (in_a) {
        m6a_in_a_gr$transcript_position[S4Vectors::subjectHits(overlaps_a)[1]]
      } else NA_integer_

      overlaps_b <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_b_gr, type = "equal")
      in_b       <- length(overlaps_b) > 0
      prob_b     <- if (in_b) {
        m6a_in_b_gr$probability[S4Vectors::subjectHits(overlaps_b)[1]]
      } else NA_real_
      tx_pos_b   <- if (in_b) {
        m6a_in_b_gr$transcript_position[S4Vectors::subjectHits(overlaps_b)[1]]
      } else NA_integer_

      transcript_pos <- if (!is.na(tx_pos_a)) tx_pos_a else tx_pos_b

      isoform_status <- if (in_a && in_b) {
        "IN_BOTH_ISOFORMS"
      } else if (in_a && !in_b) {
        "ISOFORM_A_ONLY"
      } else if (!in_a && in_b) {
        "ISOFORM_B_ONLY"
      } else {
        stop(sprintf(
          "impossible state: site %s not found in either isoform %s/%s",
          as.character(genomic_pos), iso_a, iso_b
        ))
      }

      # ── REGULATORY: which condition detected the site ───────────────────────
      cond_status <- NA_character_
      prob_c1     <- NA_real_
      prob_c2     <- NA_real_

      if (has_condition) {
        hits  <- GenomicRanges::findOverlaps(genomic_pos, m6a_sites_gr, type = "equal")
        idx   <- S4Vectors::subjectHits(hits)
        conds <- m6a_sites_gr$condition[idx]
        probs <- m6a_sites_gr$probability[idx]

        in_c1 <- !is.na(condition_1) && condition_1 %in% conds
        in_c2 <- !is.na(condition_2) && condition_2 %in% conds
        if (in_c1) prob_c1 <- probs[match(condition_1, conds)]
        if (in_c2) prob_c2 <- probs[match(condition_2, conds)]

        cond_status <- if (in_c1 && in_c2) {
          "RETAINED"
        } else if (in_c1 && !in_c2) {
          "LOST"
        } else if (!in_c1 && in_c2) {
          "GAINED"
        } else {
          NA_character_
        }
      }

      result_list[[length(result_list) + 1]] <- data.table::data.table(
        gene_id             = gene,
        isoform_a           = iso_a,
        isoform_b           = iso_b,
        condition_1         = condition_1,
        condition_2         = condition_2,
        genomic_position    = as.character(genomic_pos),
        seqname             = as.character(GenomicRanges::seqnames(genomic_pos)),
        start               = GenomicRanges::start(genomic_pos),
        end                 = GenomicRanges::end(genomic_pos),
        transcript_position = transcript_pos,
        strand              = as.character(GenomicRanges::strand(genomic_pos)),

        # structural
        m6a_in_isoform_a    = in_a,
        m6a_in_isoform_b    = in_b,
        isoform_status      = isoform_status,
        probability_a       = prob_a,
        probability_b       = prob_b,

        # regulatory
        m6a_fate            = cond_status,
        probability_c1      = prob_c1,
        probability_c2      = prob_c2,

        fdr                 = fdr,
        dif                 = dif
      )
    }
  }

  if (length(result_list) == 0) {
    warning("No m6A sites found in isoform switches. Returning empty result.",
            call. = FALSE)
    return(data.table::data.table())
  }

  result <- data.table::rbindlist(result_list)
  data.table::setorder(result, gene_id, fdr, isoform_status)

  result
}
