#' Parse m6AnetAnalyzer Output
#'
#' Reads m6A site predictions from m6AnetAnalyzer and returns a data.table
#' with standardized columns for downstream analysis.
#'
#' @param m6anet_file Path to m6AnetAnalyzer output file (typically .csv or .tsv)
#' @param probability_threshold Minimum probability to retain m6A site (default: 0.5)
#' @param transcript_col Column name for transcript ID (default: "transcript_id")
#' @param position_col Column name for genomic position (default: "position")
#' @param prob_col Column name for prediction probability (default: "probability")
#'
#' @return data.table with columns: transcript_id, position, probability, 
#'         gene_id (if available), strand, sequence_context
#'
#' @examples
#' \dontrun{
#' m6a_sites <- parse_m6anet("m6anet_predictions.csv")
#' }
#'
#' @import data.table
#' @import dplyr
#'
#' @export
parse_m6anet <- function(m6anet_file, 
                         probability_threshold = 0.5,
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
