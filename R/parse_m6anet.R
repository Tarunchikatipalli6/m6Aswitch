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
#'         and kmer (5-mer sequence context, if present in the input file).
#'         The kmer column is used by annotate_drach() for DRACH motif
#'         annotation without requiring separate transcript sequence files.
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
  data.table::setnames(m6a_data, old = c(transcript_col, position_col, prob_col),
                       new = c("transcript_id", "position", "probability"), skip_absent = TRUE)

  # Ensure required columns
  required_cols <- c("transcript_id", "position", "probability")
  if (!all(required_cols %in% names(m6a_data))) {
    stop("m6AnetAnalyzer output must contain columns: ",
         paste(required_cols, collapse = ", "))
  }

  # Filter by probability threshold
  na_prob_rows <- m6a_data[is.na(probability), .N]
  if (na_prob_rows > 0) {
    message("Dropping ", na_prob_rows, " row(s) with NA probability.")
  }
  m6a_data <- m6a_data[!is.na(probability) & probability >= probability_threshold]

  # Coerce position to numeric
  # m6Anet reports 0-based positions — convert to 1-based for downstream use
  m6a_data[, position := as.integer(position) + 1L]

  # Keep kmer column if present (used by annotate_drach() for DRACH motif annotation)
  keep_cols <- c("transcript_id", "position", "probability")
  if ("kmer" %in% names(m6a_data)) {
    keep_cols <- c(keep_cols, "kmer")
  }
  m6a_data <- m6a_data[, ..keep_cols]

  # Sort by transcript and position
  setorder(m6a_data, transcript_id, position)

  return(m6a_data)
}
