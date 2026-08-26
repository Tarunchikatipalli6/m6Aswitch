#' Classify the Mechanism Behind Isoform-Specific m6A Sites
#'
#' @md
#'
#' For every site that appears on only one isoform of a switch pair, determines
#' whether the difference is structural (the sequence is absent from the other
#' isoform) or regulatory (the sequence is present but unmethylated).
#'
#' \itemize{
#'   \item \strong{A_ONLY_EXON_SKIPPED} — the genomic position does not exist in
#'     isoform B; the exon is absent from that isoform's structure. Structural.
#'   \item \strong{A_ONLY_UNMETHYLATED} — the position exists in isoform B but
#'     m6Anet detected no modification there. Regulatory.
#'   \item \strong{B_ONLY_EXON_INCLUDED} — the position does not exist in
#'     isoform A; isoform B includes an exon that carries the site. Structural.
#'   \item \strong{B_ONLY_NEW_METHYLATION} — the position exists in isoform A but
#'     is methylated only in isoform B. Regulatory.
#' }
#'
#' Sites present on both isoforms (\code{IN_BOTH_ISOFORMS}) are left as NA.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#'   Must contain \code{seqname}, \code{start}, \code{end}, \code{isoform_a},
#'   \code{isoform_b}, and \code{isoform_status}.
#' @param gtf_file Path to the same GTF/GFF used in
#'   \code{lift_m6a_to_genomic()}.
#'
#' @return A copy of \code{m6a_switches} with two new columns:
#'   \describe{
#'     \item{isoform_mechanism}{Character. One of \code{A_ONLY_EXON_SKIPPED},
#'       \code{A_ONLY_UNMETHYLATED}, \code{B_ONLY_EXON_INCLUDED},
#'       \code{B_ONLY_NEW_METHYLATION}, or \code{NA} for sites present on both
#'       isoforms.}
#'     \item{target_isoform_has_exon}{Logical. \code{TRUE} if the genomic
#'       coordinate overlaps an exon of the isoform that \emph{lacks} the m6A
#'       call. Note the target differs by status: isoform B is checked for
#'       \code{ISOFORM_A_ONLY} sites, isoform A for \code{ISOFORM_B_ONLY}.
#'       \code{NA} for sites on both isoforms.}
#'   }
#'
#' @details
#' ## How the classification works
#'
#' For an \strong{ISOFORM_A_ONLY} site (detected on isoform A, absent on B):
#' \enumerate{
#'   \item Retrieve all exons of isoform B from the GTF.
#'   \item Check whether \code{seqname}:\code{start}-\code{end} overlaps any.
#'   \item No overlap → the exon is absent from isoform B →
#'     \code{A_ONLY_EXON_SKIPPED}.
#'   \item Overlap → the exon is present but unmethylated →
#'     \code{A_ONLY_UNMETHYLATED}.
#' }
#'
#' For an \strong{ISOFORM_B_ONLY} site the same test is applied to isoform A,
#' giving \code{B_ONLY_EXON_INCLUDED} or \code{B_ONLY_NEW_METHYLATION}.
#'
#' ## Biological interpretation
#'
#' | Classification | Type | Meaning |
#' |---|---|---|
#' | A_ONLY_EXON_SKIPPED | Structural | exon absent from isoform B; m6A cannot exist there |
#' | A_ONLY_UNMETHYLATED | Regulatory | same exon present, not methylated on isoform B |
#' | B_ONLY_EXON_INCLUDED | Structural | isoform B includes an exon carrying the site |
#' | B_ONLY_NEW_METHYLATION | Regulatory | exon present on A, methylated only on B |
#'
#' The regulatory classes cannot be explained by alternative splicing and are
#' generally the more informative.
#'
#' ## Note on interpretation
#'
#' This function operates on \code{isoform_status}, which describes isoform
#' membership rather than condition. "Unmethylated on isoform B" means m6Anet
#' made no call at that position on that transcript — which may reflect genuine
#' absence of modification, or insufficient read coverage. The function does not
#' distinguish these.
#'
#' @examples
#' \dontrun{
#' res <- annotate_m6a_switches_genomic(genomic_sites, iso_switches)
#' res <- annotate_drach(res, m6a_cond_a, m6a_cond_b)
#' res <- classify_lost_mechanism(res, "genome.gtf")
#'
#' # Regulatory changes (same sequence, different methylation)
#' res[isoform_mechanism %in% c("A_ONLY_UNMETHYLATED", "B_ONLY_NEW_METHYLATION")]
#'
#' # Structural changes (splicing)
#' res[isoform_mechanism %in% c("A_ONLY_EXON_SKIPPED", "B_ONLY_EXON_INCLUDED")]
#'
#' # Cross-tabulate against isoform status
#' table(res$isoform_status, res$isoform_mechanism, useNA = "ifany")
#' }
#'
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors queryHits
#' @import data.table
#' @import methods
#'
#' @export
classify_lost_mechanism <- function(m6a_switches, gtf_file) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!data.table::is.data.table(m6a_switches)) {
    stop("m6a_switches must be a data.table from annotate_m6a_switches_genomic()")
  }

  required_cols <- c("seqname", "start", "end", "isoform_a", "isoform_b",
                     "isoform_status")
  missing_cols  <- setdiff(required_cols, names(m6a_switches))
  if (length(missing_cols) > 0) {
    stop(
      "m6a_switches is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nThis function requires output from annotate_m6a_switches_genomic().",
      if ("m6a_fate" %in% names(m6a_switches) &&
          !"isoform_status" %in% names(m6a_switches))
        "\nFound 'm6a_fate' but not 'isoform_status' - regenerate with the current version."
      else ""
    )
  }

  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }

  # Work on a copy so the caller's table is not modified by reference
  m6a_switches <- data.table::copy(m6a_switches)

  if (nrow(m6a_switches) == 0) {
    m6a_switches[, isoform_mechanism       := character(0)]
    m6a_switches[, target_isoform_has_exon := logical(0)]
    return(m6a_switches)
  }

  # ── Build exon database ─────────────────────────────────────────────────────
  message("Building transcript exon database from GTF...")
  txdb      <- txdbmaker::makeTxDbFromGFF(gtf_file)
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  message("Exon database ready.")

  # ── Initialise output columns ───────────────────────────────────────────────
  m6a_switches[, isoform_mechanism       := NA_character_]
  m6a_switches[, target_isoform_has_exon := NA]

  # ── Classify sites present on only one isoform ──────────────────────────────
  classify_idx <- which(m6a_switches$isoform_status %in%
                        c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY"))

  if (length(classify_idx) > 0) {

    # NOTE: row_id is assigned AFTER subsetting. Inside DT[i, j], `.I` returns
    # positions within the subset, not the parent table - using it here would
    # write results to the wrong rows.
    classify_dt <- m6a_switches[classify_idx, .(
      seqname,
      start,
      end,
      isoform_status,
      target_isoform = data.table::fifelse(isoform_status == "ISOFORM_A_ONLY",
                                           isoform_b, isoform_a)
    )]
    classify_dt[, row_id := classify_idx]

    unique_targets  <- unique(classify_dt$target_isoform)
    isoform_missing <- setdiff(unique_targets, names(all_exons))
    present_targets <- intersect(unique_targets, names(all_exons))

    classify_dt[, has_exon := NA]

    for (iso in present_targets) {
      idx <- which(classify_dt$target_isoform == iso)
      site_gr <- GenomicRanges::GRanges(
        seqnames = classify_dt$seqname[idx],
        ranges   = IRanges::IRanges(start = classify_dt$start[idx],
                                    end   = classify_dt$end[idx])
      )
      hits <- GenomicRanges::findOverlaps(site_gr, all_exons[[iso]], type = "any")
      has_exon <- rep(FALSE, length(idx))
      if (length(hits) > 0) {
        has_exon[unique(S4Vectors::queryHits(hits))] <- TRUE
      }
      classify_dt$has_exon[idx] <- has_exon
    }

    if (length(isoform_missing) > 0) {
      warning(sprintf(
        "Could not find %d target isoform(s) in GTF while classifying mechanisms (%s).",
        length(isoform_missing),
        paste(utils::head(isoform_missing, 5), collapse = ", ")
      ), call. = FALSE)
    }

    m6a_switches[classify_dt$row_id, target_isoform_has_exon := classify_dt$has_exon]

    sel_a <- classify_dt$isoform_status == "ISOFORM_A_ONLY" &
             !is.na(classify_dt$has_exon)
    if (any(sel_a)) {
      m6a_switches[
        classify_dt$row_id[sel_a],
        isoform_mechanism := ifelse(
          classify_dt$has_exon[sel_a],
          "A_ONLY_UNMETHYLATED",
          "A_ONLY_EXON_SKIPPED"
        )
      ]
    }

    sel_b <- classify_dt$isoform_status == "ISOFORM_B_ONLY" &
             !is.na(classify_dt$has_exon)
    if (any(sel_b)) {
      m6a_switches[
        classify_dt$row_id[sel_b],
        isoform_mechanism := ifelse(
          classify_dt$has_exon[sel_b],
          "B_ONLY_NEW_METHYLATION",
          "B_ONLY_EXON_INCLUDED"
        )
      ]
    }
  }

  # ── Summary ─────────────────────────────────────────────────────────────────
  n_a_skip  <- sum(m6a_switches$isoform_mechanism == "A_ONLY_EXON_SKIPPED",    na.rm = TRUE)
  n_a_unmet <- sum(m6a_switches$isoform_mechanism == "A_ONLY_UNMETHYLATED",    na.rm = TRUE)
  n_b_incl  <- sum(m6a_switches$isoform_mechanism == "B_ONLY_EXON_INCLUDED",   na.rm = TRUE)
  n_b_new   <- sum(m6a_switches$isoform_mechanism == "B_ONLY_NEW_METHYLATION", na.rm = TRUE)

  message(sprintf(
    paste0("Classification complete:\n",
           "  Structural (splicing explains the difference):\n",
           "    A_ONLY_EXON_SKIPPED    : %d\n",
           "    B_ONLY_EXON_INCLUDED   : %d\n",
           "  Regulatory (position present in both isoforms):\n",
           "    A_ONLY_UNMETHYLATED    : %d\n",
           "    B_ONLY_NEW_METHYLATION : %d"),
    n_a_skip, n_b_incl, n_a_unmet, n_b_new
  ))

  m6a_switches
}
