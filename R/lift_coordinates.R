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
#' @examples
#' \dontrun{
#' m6a_sites <- parse_m6anet("m6anet_predictions.csv")
#' genomic_sites <- lift_m6a_to_genomic(m6a_sites, "genome.gtf")
#' }
#'
#' @importFrom GenomicFeatures makeTxDbFromGFF transcripts mapFromTranscripts
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
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)

  # Get all transcripts as GRanges with tx_name metadata
  all_tx <- GenomicFeatures::transcripts(txdb, columns = "tx_name")

  other_cols <- setdiff(names(m6a_sites), c("transcript_id", "position", "probability"))

  result_list <- list()

  for (i in seq_len(nrow(m6a_sites))) {
    tx_id  <- m6a_sites[i, transcript_id]
    tx_pos <- as.integer(m6a_sites[i, position])

    tx_range <- all_tx[all_tx$tx_name == tx_id]

    if (length(tx_range) == 0) {
      warning(sprintf("Transcript %s not found in GTF", tx_id))
      next
    }

    # Build a single-base query in transcript space
    tx_query <- GenomicRanges::GRanges(
      seqnames = tx_id,
      ranges   = IRanges::IRanges(start = tx_pos, end = tx_pos)
    )

    genomic_ranges <- tryCatch(
      GenomicFeatures::mapFromTranscripts(tx_query, tx_range),
      error = function(e) {
        warning(sprintf(
          "Could not map position %d in transcript %s: %s",
          tx_pos, tx_id, conditionMessage(e)
        ))
        NULL
      }
    )

    if (is.null(genomic_ranges) || length(genomic_ranges) == 0) {
      warning(sprintf(
        "Could not map position %d in transcript %s to genomic coords", tx_pos, tx_id
      ))
      next
    }

    genomic_ranges$transcript_id       <- tx_id
    genomic_ranges$transcript_position <- tx_pos
    genomic_ranges$probability         <- m6a_sites[i, probability]

    for (col in other_cols) {
      val <- m6a_sites[i, get(col)]
      if (!is.na(val)) {
        genomic_ranges[[col]] <- val
      }
    }

    result_list[[i]] <- genomic_ranges
  }

  # Remove NULL entries (from skipped transcripts)
  result_list <- Filter(Negate(is.null), result_list)

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
#' Equivalent to \code{annotate_m6a_switches()} but operates on genomic
#' coordinates (GRanges) produced by \code{lift_m6a_to_genomic()}. Because
#' m6A sites are compared at the same genomic locus across isoforms, the
#' LOST/GAINED/RETAINED calls are scientifically valid even when isoforms
#' differ in length or exon composition.
#'
#' @param m6a_sites_gr GRanges object from \code{lift_m6a_to_genomic()}.
#'   Must have metadata column \code{transcript_id} and \code{probability}.
#' @param iso_switches data.table from \code{parse_isoform_switch()}.
#' @param sequences data.table with columns \code{isoform_id} and \code{sequence}
#'   (used for DRACH annotation in downstream steps; validated here for
#'   early error detection).
#'
#' @return data.table with columns:
#'   \describe{
#'     \item{gene_id}{Gene identifier}
#'     \item{isoform_a, isoform_b}{Isoform identifiers}
#'     \item{genomic_position}{Character representation of the genomic range}
#'     \item{seqname}{Chromosome/scaffold name}
#'     \item{start, end}{Genomic coordinates (1-based)}
#'     \item{strand}{Strand of the m6A site}
#'     \item{m6a_in_isoform_a, m6a_in_isoform_b}{Logical: site detected in each isoform}
#'     \item{m6a_fate}{LOST, GAINED, or RETAINED}
#'     \item{probability_a, probability_b}{m6A probabilities per isoform}
#'     \item{fdr}{Isoform switch FDR}
#'   }
#'
#' @examples
#' \dontrun{
#' genomic_sites  <- lift_m6a_to_genomic(m6a_sites, "genome.gtf")
#' iso_switches   <- parse_isoform_switch("isoform_switches.txt")
#' iso_sequences  <- data.table::fread("isoform_sequences.csv")
#' m6a_annotated  <- annotate_m6a_switches_genomic(genomic_sites, iso_switches, iso_sequences)
#' }
#'
#' @importFrom GenomicRanges findOverlaps reduce seqnames start end strand
#' @importFrom S4Vectors subjectHits
#' @import data.table
#' @import methods
#'
#' @export
annotate_m6a_switches_genomic <- function(m6a_sites_gr, iso_switches, sequences) {

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
  if (!data.table::is.data.table(sequences)) {
    stop("sequences must be a data.table with columns: isoform_id, sequence")
  }

  if (nrow(iso_switches) == 0) {
    warning("iso_switches is empty. Returning empty result.")
    return(data.table::data.table())
  }

  result_list <- list()

  for (i in seq_len(nrow(iso_switches))) {
    switch_row <- iso_switches[i]
    gene   <- switch_row$gene_id
    iso_a  <- switch_row$isoform_a
    iso_b  <- switch_row$isoform_b
    fdr    <- switch_row$fdr

    # Subset GRanges by transcript_id metadata
    m6a_in_a_gr <- m6a_sites_gr[m6a_sites_gr$transcript_id == iso_a]
    m6a_in_b_gr <- m6a_sites_gr[m6a_sites_gr$transcript_id == iso_b]

    if (length(m6a_in_a_gr) == 0 && length(m6a_in_b_gr) == 0) {
      next
    }

    all_ranges <- c(m6a_in_a_gr, m6a_in_b_gr)
    merged     <- GenomicRanges::reduce(all_ranges, drop.empty.ranges = TRUE)

    for (j in seq_len(length(merged))) {
      genomic_pos <- merged[j]

      overlaps_a <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_a_gr, type = "equal")
      in_a       <- length(overlaps_a) > 0
      prob_a     <- if (in_a) {
        m6a_in_a_gr$probability[S4Vectors::subjectHits(overlaps_a)[1]]
      } else {
        NA_real_
      }

      overlaps_b <- GenomicRanges::findOverlaps(genomic_pos, m6a_in_b_gr, type = "equal")
      in_b       <- length(overlaps_b) > 0
      prob_b     <- if (in_b) {
        m6a_in_b_gr$probability[S4Vectors::subjectHits(overlaps_b)[1]]
      } else {
        NA_real_
      }

      fate <- if (in_a && !in_b) {
        "LOST"
      } else if (!in_a && in_b) {
        "GAINED"
      } else {
        "RETAINED"
      }

      result_list[[length(result_list) + 1]] <- data.table::data.table(
        gene_id            = gene,
        isoform_a          = iso_a,
        isoform_b          = iso_b,
        genomic_position   = as.character(genomic_pos),
        seqname            = as.character(GenomicRanges::seqnames(genomic_pos)),
        start              = GenomicRanges::start(genomic_pos),
        end                = GenomicRanges::end(genomic_pos),
        strand             = as.character(GenomicRanges::strand(genomic_pos)),
        m6a_in_isoform_a   = in_a,
        m6a_in_isoform_b   = in_b,
        m6a_fate           = fate,
        probability_a      = prob_a,
        probability_b      = prob_b,
        fdr                = fdr
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
