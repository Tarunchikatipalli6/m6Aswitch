#' Parse IsoformSwitchAnalyzeR Output
#'
#' Reads isoform switch results from IsoformSwitchAnalyzeR and returns a data.table
#' with standardized columns for integration with m6A analysis.
#'
#' @param iso_switch_file Path to IsoformSwitchAnalyzeR output (typically isoform_switches.txt)
#' @param fdr_threshold FDR threshold for significant switches (default: 0.05)
#' @param gene_col Column name for gene ID (default: "geneID")
#' @param iso_a_col Column name for isoform A (default: "isoformID_A")
#' @param iso_b_col Column name for isoform B (default: "isoformID_B")
#' @param fdr_col Column name for FDR (default: "isoform_switch_q_value")
#' @param cond_a_col Column name for condition A (default: "condition_1")
#' @param cond_b_col Column name for condition B (default: "condition_2")
#' @param dif_col Column name for dIF direction indicator (default: "dIF")
#'
#' @return data.table with columns:
#'   - gene_id, isoform_a, isoform_b, fdr
#'   - condition_1, condition_2 (condition labels for the comparison)
#'   - dif (direction and magnitude of isoform switch)
#'
#' @details
#' IsoformSwitchAnalyzeR performs pairwise comparisons between conditions.
#' Each row represents one comparison between exactly two conditions.
#' The condition_1 and condition_2 columns specify which conditions are being compared.
#'
#' @examples
#' \dontrun{
#' switches <- parse_isoform_switch(
#'   "isoform_switches.txt",
#'   cond_a_col = "condition_1",
#'   cond_b_col = "condition_2"
#' )
#' }
#'
#' @import data.table
#'
#' @export
parse_isoform_switch <- function(iso_switch_file,
                                fdr_threshold = 0.05,
                                gene_col = "geneID",
                                iso_a_col = "isoformID_A",
                                iso_b_col = "isoformID_B",
                                fdr_col = "isoform_switch_q_value",
                                cond_a_col = "condition_1",
                                cond_b_col = "condition_2",
                                dif_col = "dIF") {
  
  # Read file
  iso_data <- data.table::fread(iso_switch_file)
  
  # Standardize core column names
  data.table::setnames(iso_data,
                       old = c(gene_col, iso_a_col, iso_b_col, fdr_col),
                       new = c("gene_id", "isoform_a", "isoform_b", "fdr"),
                       skip_absent = TRUE)
  
  # Standardize condition and dIF columns
  data.table::setnames(iso_data,
                       old = c(cond_a_col, cond_b_col, dif_col),
                       new = c("condition_1", "condition_2", "dif"),
                       skip_absent = TRUE)
  
  # Filter by FDR
  iso_data <- iso_data[fdr <= fdr_threshold]
  
  # Ensure required columns
  required_cols <- c("gene_id", "isoform_a", "isoform_b", "fdr")
  if (!all(required_cols %in% names(iso_data))) {
    stop("IsoformSwitchAnalyzeR output must contain columns: ",
         paste(required_cols, collapse = ", "))
  }
  
  # Ensure condition columns exist (add NA if missing)
  if (!"condition_1" %in% names(iso_data)) {
    iso_data[, condition_1 := NA_character_]
  }
  if (!"condition_2" %in% names(iso_data)) {
    iso_data[, condition_2 := NA_character_]
  }
  if (!"dif" %in% names(iso_data)) {
    iso_data[, dif := NA_real_]
  }
  
  # Remove duplicates (keep most significant by FDR)
  iso_data <- iso_data[order(fdr)]
  iso_data <- iso_data[!duplicated(paste(gene_id, isoform_a, isoform_b))]
  
  return(iso_data)
}
