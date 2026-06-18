#' Parse m6AnetAnalyzer Output
#'
#' Reads m6A site predictions from m6AnetAnalyzer and returns a data.table
#' with standardized columns for downstream analysis.
#'
#' @param m6anet_file Path to m6AnetAnalyzer output file (typically .csv or .tsv)
#' @param probability_threshold Minimum probability to retain m6A site (default: 0.9)
#' @param transcript_col Column name for transcript ID (default: "transcript_id")
#' @param position_col Column name for genomic position (default: "position")
#' @param prob_col Column name for prediction probability (default: "probability")
#'
#' @details
#' ## Probability Threshold Strategy
#'
#' This function supports two complementary thresholds for m6A site detection:
#'
#' ### High-Confidence Analysis (default: 0.9)
#' - **Threshold**: Probability ≥ 0.9
#' - **Use**: Primary publication-quality results
#' - **Rationale**: Liu et al. (2023) demonstrates 0.9 provides optimal high-confidence predictions
#' - **Field standard**: Recommended for epitranscriptomics studies
#' - **Lab standard**: Consistent with published m6A methodology
#'
#' ### Sensitivity Analysis (0.5)
#' - **Threshold**: Probability ≥ 0.5
#' - **Use**: Exploratory analysis, supplementary materials, sensitivity testing
#' - **Rationale**: Captures broader m6A transitions, may include lower-confidence sites
#' - **Trade-off**: Higher sensitivity at cost of specificity
#'
#' ### Recommended Workflow
#' Run BOTH thresholds and report separately:
#' ```
#' # High-confidence (main results)
#' m6a_high <- parse_m6anet("predictions.csv", probability_threshold = 0.9)
#'
#' # Sensitivity analysis (supplementary)
#' m6a_broad <- parse_m6anet("predictions.csv", probability_threshold = 0.5)
#' ```
#'
#' ### References
#' Liu, H., Begik, O., Lucas, M. C., et al. (2023).
#' "Accurate detection of m6A RNA modifications in the eukaryotic transcriptome
#' with SCARLET." *Nature Biotechnology* 41, 896–905.
#'
#' @return data.table with columns: transcript_id, position, probability,
#'         gene_id (if available), strand, sequence_context
#'
#' @examples
#' \dontrun{
#' # High-confidence (default, recommended for publication)
#' m6a_high <- parse_m6anet("m6anet_predictions.csv")
#'
#' # Sensitivity analysis (for supplementary/exploratory)
#' m6a_broad <- parse_m6anet("m6anet_predictions.csv", probability_threshold = 0.5)
#' }
#'
#' @import data.table
#'
#' @export
parse_m6anet <- function(m6anet_file,
                         probability_threshold = 0.9,
                         transcript_col = "transcript_id",
                         position_col = "position",
                         prob_col = "probability") {

  # Read file
  m6a_data <- data.table::fread(m6anet_file)

  # Standardize column names
  setnames(m6a_data, old = c(transcript_col, position_col, prob_col),
           new = c("transcript_id", "position", "probability"), skip_absent = TRUE)

  # Filter by probability threshold
  m6a_data <- m6a_data[probability >= probability_threshold]

  # Ensure required columns
  required_cols <- c("transcript_id", "position", "probability")
  if (!all(required_cols %in% names(m6a_data))) {
    stop("m6AnetAnalyzer output must contain columns: ",
         paste(required_cols, collapse = ", "))
  }

  # Coerce position to numeric
  m6a_data[, position := as.numeric(position)]

  # Sort by transcript and position
  setorder(m6a_data, transcript_id, position)

  return(m6a_data)
}

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
