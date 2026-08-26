#' Parse IsoformSwitchAnalyzeR Output
#'
#' Reads isoform switch pairs and returns a data.table with standardised
#' columns for integration with m6A analysis.
#'
#' @param iso_switch_file Path to a switch-pair table. Typically produced from
#'   \code{isoformSwitchTestDEXSeq()} output by pairing the up- and
#'   down-regulated isoform of each gene.
#' @param fdr_threshold FDR cutoff for significant switches (default 0.05).
#' @param gene_col Column name for gene ID (default "geneID").
#' @param iso_a_col Column name for isoform A (default "isoformID_A").
#' @param iso_b_col Column name for isoform B (default "isoformID_B").
#' @param fdr_col Column name for FDR (default "isoform_switch_q_value").
#' @param cond_a_col Column name for condition 1 (default "condition_1").
#' @param cond_b_col Column name for condition 2 (default "condition_2").
#' @param dif_col Column name for dIF (default "dIF").
#' @param verbose Logical. Print a summary of the parsed comparison, including
#'   the direction of the A/B assignment (default TRUE).
#'
#' @return data.table with columns \code{gene_id}, \code{isoform_a},
#'   \code{isoform_b}, \code{fdr}, \code{condition_1}, \code{condition_2},
#'   \code{dif}. One row per switch pair, deduplicated on the gene-isoform
#'   triple, keeping the most significant.
#'
#' @section What isoform A and B mean:
#' IsoformSwitchAnalyzeR defines \code{dIF = IF(condition_2) - IF(condition_1)}.
#' The switch-pair table is normally built as:
#' \preformatted{
#' up   <- sub[dIF > 0]   # gaining usage in condition_2
#' down <- sub[dIF < 0]   # losing usage in condition_2
#' isoformID_A = down$isoform_id[1]
#' isoformID_B = up$isoform_id[1]
#' }
#' So \strong{isoform A is favoured in condition_1} and \strong{isoform B is
#' favoured in condition_2}.
#'
#' Condition order is set by factor level in IsoformSwitchAnalyzeR, which may
#' not be the order you expect - \code{"IDH_MUT"} sorts before \code{"IDH_WT"},
#' for instance, making the mutant condition_1. Check the printed summary, or:
#' \preformatted{
#' head(iso_switches[, .(condition_1, condition_2, dif)])
#' }
#'
#' To control it, set the factor level explicitly in the design matrix before
#' running IsoformSwitchAnalyzeR:
#' \preformatted{
#' design$condition <- factor(design$condition, levels = c("WT", "MUT"))
#' }
#'
#' @examples
#' \dontrun{
#' switches <- parse_isoform_switch("switch_pairs.txt")
#'
#' # Custom column names
#' switches <- parse_isoform_switch(
#'   "switches.txt",
#'   gene_col  = "gene_name",
#'   iso_a_col = "iso_down",
#'   iso_b_col = "iso_up"
#' )
#' }
#'
#' @import data.table
#'
#' @export
parse_isoform_switch <- function(iso_switch_file,
                                 fdr_threshold = 0.05,
                                 gene_col   = "geneID",
                                 iso_a_col  = "isoformID_A",
                                 iso_b_col  = "isoformID_B",
                                 fdr_col    = "isoform_switch_q_value",
                                 cond_a_col = "condition_1",
                                 cond_b_col = "condition_2",
                                 dif_col    = "dIF",
                                 verbose    = TRUE) {

  if (!file.exists(iso_switch_file)) {
    stop("Isoform switch file not found: ", iso_switch_file)
  }
  if (!is.numeric(fdr_threshold) || fdr_threshold < 0 || fdr_threshold > 1) {
    stop("fdr_threshold must be between 0 and 1")
  }

  iso_data <- data.table::fread(iso_switch_file)

  if (nrow(iso_data) == 0) {
    stop("Isoform switch file is empty: ", iso_switch_file)
  }

  n_input <- nrow(iso_data)

  # Standardise core column names
  data.table::setnames(
    iso_data,
    old = c(gene_col, iso_a_col, iso_b_col, fdr_col),
    new = c("gene_id", "isoform_a", "isoform_b", "fdr"),
    skip_absent = TRUE
  )

  # Standardise condition and dIF columns
  data.table::setnames(
    iso_data,
    old = c(cond_a_col, cond_b_col, dif_col),
    new = c("condition_1", "condition_2", "dif"),
    skip_absent = TRUE
  )

  required_cols <- c("gene_id", "isoform_a", "isoform_b", "fdr")
  missing_cols  <- setdiff(required_cols, names(iso_data))
  if (length(missing_cols) > 0) {
    stop("Isoform switch file is missing columns: ",
         paste(missing_cols, collapse = ", "),
         "\n  Found: ", paste(names(iso_data), collapse = ", "),
         "\n  Set gene_col / iso_a_col / iso_b_col / fdr_col to match your file.")
  }

  # Filter by FDR
  na_fdr_rows <- iso_data[is.na(fdr), .N]
  if (na_fdr_rows > 0) {
    message("Dropping ", na_fdr_rows, " row(s) with NA FDR.")
  }
  iso_data <- iso_data[!is.na(fdr) & fdr <= fdr_threshold]

  if (nrow(iso_data) == 0) {
    warning("No switches passed fdr_threshold = ", fdr_threshold,
            " (input had ", n_input, " rows).", call. = FALSE)
  }

  # Fill missing optional columns
  if (!"condition_1" %in% names(iso_data)) {
    iso_data[, condition_1 := NA_character_]
  }
  if (!"condition_2" %in% names(iso_data)) {
    iso_data[, condition_2 := NA_character_]
  }
  if (!"dif" %in% names(iso_data)) {
    iso_data[, dif := NA_real_]
  }

  # Deduplicate on the gene-isoform triple, keeping the most significant
  n_before <- nrow(iso_data)
  data.table::setorder(iso_data, fdr)
  iso_data <- unique(iso_data, by = c("gene_id", "isoform_a", "isoform_b"))
  n_dropped <- n_before - nrow(iso_data)
  if (n_dropped > 0) {
    message("Removed ", n_dropped, " duplicate switch pair(s), ",
            "keeping the most significant.")
  }

  # ── Summary, including the A/B direction ────────────────────────────────────
  if (isTRUE(verbose) && nrow(iso_data) > 0) {
    message(sprintf("Parsed %s: %d switch pair(s) at FDR <= %.3g",
                    basename(iso_switch_file), nrow(iso_data), fdr_threshold))

    c1 <- unique(iso_data$condition_1)
    c2 <- unique(iso_data$condition_2)

    if (length(c1) == 1 && length(c2) == 1 && !is.na(c1) && !is.na(c2)) {
      message(sprintf("  Comparison: %s (condition_1) vs %s (condition_2)", c1, c2))
      message(sprintf("  isoform_a is favoured in %s", c1))
      message(sprintf("  isoform_b is favoured in %s", c2))
    } else if (length(c1) > 1 || length(c2) > 1) {
      message("  Multiple comparisons present: ",
              paste(unique(paste(iso_data$condition_1, "vs", iso_data$condition_2)),
                    collapse = "; "))
    } else {
      message("  No condition labels found - A/B direction cannot be reported.")
    }

    # Sanity check: dIF should be positive if B was taken as the up isoform
    if ("dif" %in% names(iso_data) && any(!is.na(iso_data$dif))) {
      n_pos <- sum(iso_data$dif > 0, na.rm = TRUE)
      n_neg <- sum(iso_data$dif < 0, na.rm = TRUE)
      if (n_neg > n_pos) {
        warning("Most dIF values are negative (", n_neg, " negative vs ",
                n_pos, " positive). The dIF column normally records the ",
                "up-regulated isoform (B), so this suggests the A/B ",
                "assignment may be reversed relative to the usual convention. ",
                "Verify against IF1/IF2 in your isoformFeatures table.",
                call. = FALSE)
      }
    }

    gene_counts <- iso_data[, .N, by = gene_id]
    n_multi <- sum(gene_counts$N > 1)
    if (n_multi > 0) {
      message("  ", n_multi, " gene(s) have more than one switch pair.")
    }
  }

  iso_data
}
