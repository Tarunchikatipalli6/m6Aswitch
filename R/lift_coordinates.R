#' Lift m6A Transcript Coordinates to Genomic Coordinates
#'
#' Converts transcript-level m6A coordinates (as reported by m6AnetAnalyzer)
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
#' The GTF file serves as a map that specifies which exons belong to each transcript
#' and their exact positions on the chromosome. This function uses that map to
#' convert from transcript numbering to chromosome numbering.
#'
#' The key insight: mapFromTranscripts() needs the actual exon structure to correctly
#' skip over introns. Using transcript boundaries alone would place position 150 in
#' the intron (wrong). Using exon-by-exon structure allows correct conversion:
#'
#' Example:
#'   Isoform A exon structure (from GTF):
#'     Exon1: chr10:1000-1100 (positions 1-100 in transcript)
#'     Exon2: chr10:5000-5200 (positions 101-300 in transcript)
#'     [Intron: chr10:1100-5000 — skipped]
#'
#'   Input: position 150 in transcript ENST001
#'   Lookup: Position 150 is 50bp into Exon2
#'   Calculation: Exon2 starts at chr10:5000, so position 150 = chr10:5050
#'   Output: GRanges with seqname=chr10, start=5050, end=5050
#'
#' If we only used transcript boundaries (start=1000, end=5200):
#'   Naive calculation: 1000 + 150 = chr10:1150 ✗ WRONG (inside intron)
#'
#' With exon structure:
#'   Correct calculation: chr10:5050 ✓ RIGHT (in Exon2)
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

  # Get exon structure for each transcript (needed for correct coordinate conversion)
  # exonsBy() returns the actual exons that make up each transcript
  # This is critical because mapFromTranscripts() needs to know where introns are
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)

  other_cols <- setdiff(names(m6a_sites), c("transcript_id", "position", "probability"))

  result_list <- list()
  missing_tx <- data.table::data.table(transcript_id = character(), n = integer())
  failed_map <- data.table::data.table(transcript_id = character(), n = integer())

  sites_by_tx <- split(m6a_sites, by = "transcript_id", keep.by = FALSE, drop = TRUE)

  for (tx_id in names(sites_by_tx)) {
    tx_sites <- sites_by_tx[[tx_id]]
    tx_range <- all_exons[tx_id]

    if (is.null(tx_range) || length(tx_range) == 0) {
      missing_tx <- rbind(
        missing_tx,
        data.table::data.table(transcript_id = tx_id, n = nrow(tx_sites))
      )
      next
    }

    tx_pos <- as.integer(tx_sites$position)
    valid_pos <- !is.na(tx_pos)
    if (!all(valid_pos)) {
      failed_map <- rbind(
        failed_map,
        data.table::data.table(transcript_id = tx_id, n = sum(!valid_pos))
      )
    }
    if (!any(valid_pos)) next

    tx_sites_valid <- tx_sites[valid_pos]
    tx_pos_valid <- tx_pos[valid_pos]

    tx_query <- GenomicRanges::GRanges(
      seqnames = rep(tx_id, length(tx_pos_valid)),
      ranges = IRanges::IRanges(start = tx_pos_valid, end = tx_pos_valid)
    )

    genomic_ranges <- tryCatch(
      GenomicFeatures::mapFromTranscripts(tx_query, tx_range),
      error = function(e) {
        warning("mapFromTranscripts failed for ", tx_id, ": ", conditionMessage(e),
            call. = FALSE)
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
        data.table::data.table(transcript_id = tx_id, n = nrow(tx_sites_valid) - nrow(mapped_sites))
      )
    }

    genomic_ranges$transcript_id <- mapped_sites$transcript_id
    genomic_ranges$transcript_position <- as.integer(mapped_sites$position)
    genomic_ranges$probability <- mapped_sites$probability

    for (col in other_cols) {
      genomic_ranges[[col]] <- mapped_sites[[col]]
    }

    result_list[[length(result_list) + 1L]] <- genomic_ranges
  }

  # Remove NULL entries (from skipped transcripts)
  result_list <- Filter(Negate(is.null), result_list)

  if (nrow(missing_tx) > 0) {
    miss_sum <- missing_tx[, .(n = sum(n)), by = transcript_id]
    warning(
      sprintf(
        "Skipped %d site(s): %d transcript(s) not found in GTF (%s).",
        sum(miss_sum$n),
        nrow(miss_sum),
        paste(head(miss_sum$transcript_id, 5), collapse = ", ")
      )
    )
  }
  if (nrow(failed_map) > 0) {
    fail_sum <- failed_map[, .(n = sum(n)), by = transcript_id]
    warning(
      sprintf(
        "Could not map %d site(s) across %d transcript(s).",
        sum(fail_sum$n),
        nrow(fail_sum)
      )
    )
  }

  if (length(result_list) == 0) {
    stop(paste(
      "No m6A sites could be mapped to genomic coordinates.",
      "Check GTF file and transcript IDs."
    ))
  }

  result <- do.call(c, result_list)
  return(result)
}


#' Annotate Isoform Switches Using Genomic m6A Coordinates
#'
#' RECOMMENDED FOR PUBLICATION USE.
#'
#' Equivalent to \code{annotate_m6a_switches()} but operates on genomic
#' coordinates (GRanges) produced by \code{lift_m6a_to_genomic()}. Because
#' m6A sites are compared at the same genomic locus across isoforms, the
#' LOST/GAINED/RETAINED calls are scientifically valid even when isoforms
#' differ in length or exon composition.
#'
#' @param m6a_sites_gr GRanges object from \code{lift_m6a_to_genomic()}.
#'   Must have metadata column \code{transcript_id} and \code{probability}.
#' @param iso_switches data.table from \code{parse_isoform_switch()}.
#'   Should include condition_1, condition_2, and dif columns for full
#'   biological interpretation.
#'
#' @return data.table with columns:
#'   \describe{
#'     \item{gene_id}{Gene identifier}
#'     \item{isoform_a, isoform_b}{Isoform identifiers}
#'     \item{condition_1, condition_2}{Condition labels (pairwise comparison)}
#'     \item{genomic_position}{Character representation of the genomic range}
#'     \item{seqname}{Chromosome/scaffold name}
#'     \item{start, end}{Genomic coordinates (1-based) for comparison}
#'     \item{transcript_position}{Original transcript coordinate (for motif lookup)}
#'     \item{strand}{Strand of the m6A site}
#'     \item{m6a_in_isoform_a, m6a_in_isoform_b}{Logical: site detected in each isoform}
#'     \item{m6a_fate}{LOST, GAINED, or RETAINED}
#'     \item{probability_a, probability_b}{m6A probabilities per isoform}
#'     \item{fdr}{Isoform switch FDR}
#'     \item{dif}{Direction indicator from isoform switch}
#'   }
#'
#' @details
#' This is the publication-quality version of m6A-isoform integration because:
#'
#' 1. GENOMIC COORDINATES: Compares m6A sites at the same genomic location across
#'    isoforms, not just matching transcript position numbers.
#'
#' 2. EXON STRUCTURE AWARE: Correctly handles isoforms with different exon structures
#'    (e.g., one isoform skips an exon). Position 150 in two isoforms with different
#'    exon structures will be compared at their actual genomic locations, which may
#'    be completely different.
#'
#' 3. FALSE CALL PREVENTION: Avoids false LOST/GAINED/RETAINED calls that arise from
#'    comparing transcript position numbers across isoforms with different structures.
#'
#' 4. DRACH VALIDATION READY: Preserves transcript_position for motif lookup in
#'    downstream annotate_drach() function.
#'
#' Example of why this matters:
#'   Isoform A: Exon1 (chr10:1000-1100) + Exon2 (chr10:5000-5200)
#'   Isoform B: Exon1 (chr10:1000-1100) + Exon3 (chr10:8000-8200)  [Exon2 skipped]
#'
#'   Using transcript positions (WRONG):
#'     Position 150 in A: has m6A
#'     Position 150 in B: has m6A
#'     Tool says: RETAINED ✗ (FALSE)
#'
#'   Using genomic coordinates (CORRECT):
#'     Genomic chr10:5050 in A: has m6A
#'     Genomic chr10:5050 in B: DOESN'T EXIST (region skipped)
#'     Tool says: LOST ✓ (TRUE)
#'
#' m6A fate categories (genomically valid):
#' - LOST: m6A at this genomic location in isoform_a but not in isoform_b
#' - GAINED: m6A at this genomic location in isoform_b but not in isoform_a
#' - RETAINED: m6A at this genomic location in both isoforms
#'
#' @seealso annotate_m6a_switches() for quick exploratory analysis (not for publication)
#'
#' @examples
#' \dontrun{
#' # PUBLICATION-QUALITY WORKFLOW
#' m6a_sites <- parse_m6anet("m6anet_predictions.csv")
#' iso_switches <- parse_isoform_switch("isoform_switches.txt")
#' # Step 1: Lift to genomic coordinates
#' genomic_sites <- lift_m6a_to_genomic(m6a_sites, "genome.gtf")
#'
#' # Step 2: Compare using genomic coordinates (scientifically valid)
#' m6a_switches <- annotate_m6a_switches_genomic(
#'   genomic_sites, iso_switches
#' )
#'
#' # Step 3: Add functional annotation
#' m6a_final <- annotate_drach(m6a_switches, m6a_condition_a, m6a_condition_b)
#'
#' # Step 4: Report results
#' # Now your LOST/GAINED/RETAINED calls are biologically defensible
#' }
#'
#' @importFrom GenomicRanges findOverlaps seqnames start end strand
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
    warning("iso_switches is empty. Returning empty result.")
    return(data.table::data.table())
  }

  m6a_by_tx <- split(m6a_sites_gr, m6a_sites_gr$transcript_id)
  result_list <- list()

  for (i in seq_len(nrow(iso_switches))) {
    switch_row <- iso_switches[i]
    gene   <- switch_row$gene_id
    iso_a  <- switch_row$isoform_a
    iso_b  <- switch_row$isoform_b
    fdr    <- switch_row$fdr

    # Extract condition information (pairwise comparison)
    condition_1 <- if ("condition_1" %in% names(switch_row))
      switch_row$condition_1 else NA_character_
    condition_2 <- if ("condition_2" %in% names(switch_row))
      switch_row$condition_2 else NA_character_

    # Extract dIF (direction indicator)
    dif <- if ("dif" %in% names(switch_row))
      switch_row$dif else NA_real_

    # Subset GRanges by transcript_id metadata (genomic coordinates)
    m6a_in_a_gr <- m6a_by_tx[[iso_a]]
    m6a_in_b_gr <- m6a_by_tx[[iso_b]]
    if (is.null(m6a_in_a_gr)) m6a_in_a_gr <- GenomicRanges::GRanges()
    if (is.null(m6a_in_b_gr)) m6a_in_b_gr <- GenomicRanges::GRanges()

    if (length(m6a_in_a_gr) == 0 && length(m6a_in_b_gr) == 0) {
      next
    }

    # Keep distinct single-base genomic loci (do not merge adjacent sites)
    all_ranges <- c(m6a_in_a_gr, m6a_in_b_gr)
    merged     <- unique(all_ranges)

    # For each unique genomic position, determine fate
    for (j in seq_len(length(merged))) {
      genomic_pos <- merged[j]

      # Check if this genomic position has m6A in isoform A
      overlaps_a <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_a_gr, type = "equal")
      in_a       <- length(overlaps_a) > 0
      prob_a     <- if (in_a) {
        m6a_in_a_gr$probability[S4Vectors::subjectHits(overlaps_a)[1]]
      } else {
        NA_real_
      }

      # Get transcript position from isoform A (if available)
      tx_pos_a <- if (in_a) {
        m6a_in_a_gr$transcript_position[S4Vectors::subjectHits(overlaps_a)[1]]
      } else {
        NA_integer_
      }

      # Check if this genomic position has m6A in isoform B
      overlaps_b <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_b_gr, type = "equal")
      in_b       <- length(overlaps_b) > 0
      prob_b     <- if (in_b) {
        m6a_in_b_gr$probability[S4Vectors::subjectHits(overlaps_b)[1]]
      } else {
        NA_real_
      }

      # Get transcript position from isoform B (if available)
      tx_pos_b <- if (in_b) {
        m6a_in_b_gr$transcript_position[S4Vectors::subjectHits(overlaps_b)[1]]
      } else {
        NA_integer_
      }

      # Use the transcript position that's available (prefer non-NA)
      transcript_pos <- if (!is.na(tx_pos_a)) tx_pos_a else tx_pos_b

      # Determine fate (genomically valid)
      fate <- if (in_a && in_b) {
        "RETAINED"
      } else if (in_a && !in_b) {
        "LOST"       # m6A at this genomic locus in condition_1 but not condition_2
      } else if (!in_a && in_b) {
        "GAINED"     # m6A at this genomic locus not in condition_1 but present in condition_2
      } else {
        stop(
          sprintf(
            "impossible state: site %s not found in either isoform %s/%s",
            as.character(genomic_pos), iso_a, iso_b
          )
        )
      }

      result_list[[length(result_list) + 1]] <- data.table::data.table(
        gene_id            = gene,
        isoform_a          = iso_a,
        isoform_b          = iso_b,
        condition_1        = condition_1,
        condition_2        = condition_2,
        genomic_position   = as.character(genomic_pos),
        seqname            = as.character(GenomicRanges::seqnames(genomic_pos)),
        start              = GenomicRanges::start(genomic_pos),
        end                = GenomicRanges::end(genomic_pos),
        transcript_position = transcript_pos,
        strand             = as.character(GenomicRanges::strand(genomic_pos)),
        m6a_in_isoform_a   = in_a,
        m6a_in_isoform_b   = in_b,
        m6a_fate           = fate,
        probability_a      = prob_a,
        probability_b      = prob_b,
        fdr                = fdr,
        dif                = dif
      )
    }
  }

  if (length(result_list) == 0) {
    warning("No m6A sites found in isoform switches. Returning empty result.")
    return(data.table::data.table())
  }

  result <- data.table::rbindlist(result_list)
  data.table::setorder(result, gene_id, fdr, m6a_fate)

  return(result)
}
