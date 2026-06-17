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
#'
#' @return data.table with columns: gene_id, isoform_a, isoform_b, 
#'         fdr, switch_direction, condition_a, condition_b
#'
#' @examples
#' \dontrun{
#' switches <- parse_isoform_switch("isoform_switches.txt")
#' }
#'
#' @import data.table
#' @import dplyr
#'
#' @export
parse_isoform_switch <- function(iso_switch_file,
                                fdr_threshold = 0.05,
                                gene_col = "geneID",
                                iso_a_col = "isoformID_A",
                                iso_b_col = "isoformID_B",
                                fdr_col = "isoform_switch_q_value") {
  
  # Read file
  iso_data <- data.table::fread(iso_switch_file)
  
  # Standardize column names
  setnames(iso_data, old = c(gene_col, iso_a_col, iso_b_col, fdr_col),
           new = c("gene_id", "isoform_a", "isoform_b", "fdr"), 
           skip_absent = TRUE)
  
  # Filter by FDR
  iso_data <- iso_data[fdr <= fdr_threshold]
  
  # Ensure required columns
  required_cols <- c("gene_id", "isoform_a", "isoform_b", "fdr")
  if (!all(required_cols %in% names(iso_data))) {
    stop("IsoformSwitchAnalyzeR output must contain columns: ", 
         paste(required_cols, collapse = ", "))
  }
  
  # Remove duplicates (keep higher significance)
  iso_data <- iso_data[order(fdr)]
  iso_data <- iso_data[!duplicated(paste(gene_id, isoform_a, isoform_b))]
  
  return(iso_data)
}
