#' Classify the Mechanism Behind LOST m6A Sites
#'
#' @md
#'
#' For every LOST call in an annotated m6A switches table, determines whether the
#' site was lost because:
#' \itemize{
#'   \item \strong{LOST_EXON_SKIPPED} — the genomic region containing the m6A site
#'     does not exist in isoform B (the exon is absent from that isoform's exon
#'     structure).
#'   \item \strong{LOST_UNMETHYLATED} — the exon is present in isoform B but
#'     m6Anet did not detect m6A there (the sequence exists but the methylation
#'     machinery chose not to deposit m6A at that site).
#' }
#' GAINED sites receive an analogous classification (\code{GAINED_EXON_INCLUDED}
#' vs \code{GAINED_NEW_METHYLATION}). RETAINED sites are left unchanged.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()} (or
#'   any table that contains columns \code{seqname}, \code{start}, \code{end},
#'   \code{isoform_b} (for LOST) / \code{isoform_a} (for GAINED), and
#'   \code{m6a_fate}).
#' @param gtf_file Path to the same GTF/GFF file used in
#'   \code{lift_m6a_to_genomic()}.
#'
#' @return \code{m6a_switches} with two new columns added in-place:
#'   \describe{
#'     \item{lost_mechanism}{Character. One of \code{LOST_EXON_SKIPPED},
#'       \code{LOST_UNMETHYLATED}, \code{GAINED_EXON_INCLUDED},
#'       \code{GAINED_NEW_METHYLATION}, or \code{NA} for RETAINED sites.}
#'     \item{isoform_b_has_exon}{Logical. \code{TRUE} if the genomic coordinate of
#'       the m6A site overlaps at least one exon of the isoform that \emph{lacks}
#'       the m6A signal. \code{NA} for RETAINED sites.}
#'   }
#'
#' @details
#' ## How the classification works
#'
#' For a \strong{LOST} site (m6A detected in isoform A, absent in isoform B):
#' \enumerate{
#'   \item The function retrieves all exons of isoform B from the GTF.
#'   \item It checks whether the genomic coordinate (\code{seqname}:\code{start}-
#'     \code{end}) of the LOST site overlaps any of those exons.
#'   \item If \strong{no overlap} → the exon containing the m6A site does not
#'     exist in isoform B → \code{LOST_EXON_SKIPPED}.
#'   \item If \strong{overlap} → the exon exists in isoform B but m6Anet found no
#'     m6A there → \code{LOST_UNMETHYLATED}.
#' }
#'
#' For a \strong{GAINED} site (m6A absent in isoform A, detected in isoform B)
#' the same logic is applied to isoform A:
#' \enumerate{
#'   \item If the region is \strong{absent} from isoform A's exons →
#'     \code{GAINED_EXON_INCLUDED} (a new exon appears in isoform B that brings
#'     the m6A site along).
#'   \item If the region is \strong{present} in isoform A but unmethylated →
#'     \code{GAINED_NEW_METHYLATION} (same exon, but methylation is deposited only
#'     in isoform B).
#' }
#'
#' ## Biological interpretation
#'
#' | Classification | Meaning |
#' |---|---|
#' | LOST_EXON_SKIPPED | Structural — the exon is absent from isoform B; m6A physically cannot exist there |
#' | LOST_UNMETHYLATED | Regulatory — same exon present, but METTL3 does not methylate the site in isoform B |
#' | GAINED_EXON_INCLUDED | Structural — a new exon appears in isoform B that carries an m6A site |
#' | GAINED_NEW_METHYLATION | Regulatory — the exon exists in isoform A but becomes methylated only in isoform B |
#'
#' \code{LOST_UNMETHYLATED} and \code{GAINED_NEW_METHYLATION} are the most
#' biologically interesting because they represent genuine changes in m6A
#' deposition that cannot be explained by alternative splicing alone.
#'
#' @examples
#' \dontrun{
#' # After the standard publication workflow:
#' m6a_final <- annotate_drach(m6a_switches, m6a_condition_a, m6a_condition_b)
#'
#' # Classify why each LOST/GAINED site changed:
#' m6a_classified <- classify_lost_mechanism(m6a_final, "genome.gtf")
#'
#' # Inspect regulatory changes (same exon, different methylation)
#' m6a_classified[lost_mechanism == "LOST_UNMETHYLATED"]
#' m6a_classified[lost_mechanism == "GAINED_NEW_METHYLATION"]
#'
#' # Inspect structural changes (exon gained/lost)
#' m6a_classified[lost_mechanism == "LOST_EXON_SKIPPED"]
#' m6a_classified[lost_mechanism == "GAINED_EXON_INCLUDED"]
#' }
#'
#' @importFrom txdbmaker makeTxDbFromGFF
#' @importFrom GenomicFeatures exonsBy
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom IRanges IRanges
#' @import data.table
#' @import methods
#'
#' @export
classify_lost_mechanism <- function(m6a_switches, gtf_file) {

  # ── Input validation ────────────────────────────────────────────────────────
  if (!data.table::is.data.table(m6a_switches)) {
    stop("m6a_switches must be a data.table from annotate_m6a_switches_genomic()")
  }

  required_cols <- c("seqname", "start", "end", "isoform_a", "isoform_b", "isoform_status")
  missing_cols  <- setdiff(required_cols, names(m6a_switches))
  if (length(missing_cols) > 0) {
    stop(
      "m6a_switches is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      "\nThis function requires output from annotate_m6a_switches_genomic()."
    )
  }

  if (!file.exists(gtf_file)) {
    stop("GTF file not found: ", gtf_file)
  }

  if (nrow(m6a_switches) == 0) {
    m6a_switches[, isoform_mechanism        := character(0)]
    m6a_switches[, target_isoform_has_exon  := logical(0)]
    return(m6a_switches)
  }

  # ── Build exon database ──────────────────────────────────────────────────────
  message("Building transcript exon database from GTF...")
  txdb      <- txdbmaker::makeTxDbFromGFF(gtf_file)
  all_exons <- GenomicFeatures::exonsBy(txdb, by = "tx", use.names = TRUE)
  message("Exon database ready.")

  # ── Initialise output columns ────────────────────────────────────────────────
  m6a_switches[, isoform_mechanism       := NA_character_]
  m6a_switches[, target_isoform_has_exon := NA]

  # ── Vectorised classification by unique site/transcript target ───────────────
  classify_idx <- which(m6a_switches$isoform_status %in%
                        c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY"))

  if (length(classify_idx) > 0) {
    classify_dt <- m6a_switches[classify_idx, .(
      seqname,
      start,
      end,
      isoform_status,
      target_isoform = data.table::fifelse(isoform_status == "ISOFORM_A_ONLY",
                                           isoform_b, isoform_a)
    )]
    classify_dt[, row_id := classify_idx]

    unique_targets <- unique(classify_dt$target_isoform)
    isoform_missing <- setdiff(unique_targets, names(all_exons))

    classify_dt[, has_exon := NA]
    present_targets <- intersect(unique_targets, names(all_exons))

    for (iso in present_targets) {
      idx <- which(classify_dt$target_isoform == iso)
      site_gr <- GenomicRanges::GRanges(
        seqnames = classify_dt$seqname[idx],
        ranges = IRanges::IRanges(start = classify_dt$start[idx], end = classify_dt$end[idx])
      )
      hits <- GenomicRanges::findOverlaps(site_gr, all_exons[[iso]], type = "any")
      has_exon <- rep(FALSE, length(idx))
      if (length(hits) > 0) {
        has_exon[unique(S4Vectors::queryHits(hits))] <- TRUE
      }
      classify_dt$has_exon[idx] <- has_exon
    }

    if (length(isoform_missing) > 0) {
      warning(
        sprintf(
          "Could not find %d target isoform(s) in GTF while classifying mechanisms.",
          length(isoform_missing)
        )
      )
    }

    m6a_switches[classify_dt$row_id, target_isoform_has_exon := classify_dt$has_exon]

    sel_a <- classify_dt$isoform_status == "ISOFORM_A_ONLY" & !is.na(classify_dt$has_exon)
    m6a_switches[
      classify_dt$row_id[sel_a],
      isoform_mechanism := ifelse(
        classify_dt$has_exon[sel_a],
        "A_ONLY_UNMETHYLATED",
        "A_ONLY_EXON_SKIPPED"
      )
    ]

    sel_b <- classify_dt$isoform_status == "ISOFORM_B_ONLY" & !is.na(classify_dt$has_exon)
    m6a_switches[
      classify_dt$row_id[sel_b],
      isoform_mechanism := ifelse(
        classify_dt$has_exon[sel_b],
        "B_ONLY_NEW_METHYLATION",
        "B_ONLY_EXON_INCLUDED"
      )
    ]
  }

  # ── Summary message ─────────────────────────────────────────────────────────
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

  return(m6a_switches)
}
