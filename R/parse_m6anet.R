#' Parse m6Anet Output
#'
#' Reads m6A site predictions from m6Anet and returns a data.table with
#' standardised columns for downstream analysis.
#'
#' @param m6anet_file Path to m6Anet output (typically \code{data.site_proba.csv}).
#' @param probability_threshold Minimum probability to retain a site
#'   (default 0.9). Set to 0 to retain every position m6Anet evaluated - see
#'   the testability section below.
#' @param transcript_col Column name for transcript ID (default "transcript_id").
#' @param position_col Column name for transcript position (default "position").
#' @param prob_col Column name for modification probability (default "probability").
#' @param keep_extra Character vector of additional columns to retain. Useful
#'   for carrying \code{n_reads}, \code{mod_ratio}, or a user-added
#'   \code{condition} column through to later steps. \code{kmer} is retained
#'   automatically when present.
#'
#' @return data.table with columns \code{transcript_id}, \code{position},
#'   \code{probability}, \code{kmer} (when present), plus anything named in
#'   \code{keep_extra}. Sorted by transcript and position.
#'
#' @details
#' ## Coordinate convention
#'
#' m6Anet reports **0-based** transcript positions. This function converts them
#' to 1-based to match R and Bioconductor. Verified against transcript sequence:
#' at the raw position the base is usually G, at position + 1 it is always A,
#' and the reported \code{kmer} matches the sequence only after the shift.
#'
#' ## Probability threshold
#'
#' | Threshold | Use |
#' |---|---|
#' | 0.9 | Primary, publication-quality |
#' | 0.8 | Moderate |
#' | 0.5 | Sensitivity analysis, supplementary |
#' | 0.0 | Every evaluated position - for building a testability filter |
#'
#' Running both a strict and a permissive threshold and reporting them
#' separately is good practice.
#'
#' ## Testability filtering
#'
#' m6Anet only reports a position when read coverage is sufficient. A position
#' absent from one condition's output may be unmethylated, or simply untested.
#' Treating the two as equivalent produces spurious gains when one library is
#' shallower.
#'
#' To distinguish them, parse both conditions at \code{probability_threshold = 0}
#' and intersect the positions:
#' \preformatted{
#' a_all <- parse_m6anet("cond_a.csv", probability_threshold = 0)
#' b_all <- parse_m6anet("cond_b.csv", probability_threshold = 0)
#'
#' testable <- merge(unique(a_all[, .(transcript_id, position)]),
#'                   unique(b_all[, .(transcript_id, position)]),
#'                   by = c("transcript_id", "position"))
#'
#' a <- merge(parse_m6anet("cond_a.csv", 0.9), testable,
#'            by = c("transcript_id", "position"))
#' b <- merge(parse_m6anet("cond_b.csv", 0.9), testable,
#'            by = c("transcript_id", "position"))
#' }
#'
#' ## Carrying condition labels
#'
#' To enable condition-level comparison downstream, tag each table before
#' combining:
#' \preformatted{
#' a[, condition := "WT"]
#' b[, condition := "MUT"]
#' combined <- rbind(a, b)   # do NOT deduplicate across conditions
#' }
#'
#' ## Reference
#'
#' Liu, H., Begik, O., Lucas, M. C., et al. (2023). Accurate detection of m6A
#' RNA modifications in the eukaryotic transcriptome with SCARLET.
#' *Nature Biotechnology* 41, 896-905.
#'
#' @examples
#' \dontrun{
#' # Default column names
#' m6a <- parse_m6anet("m6anet_predictions.csv")
#'
#' # m6Anet's own output naming
#' m6a <- parse_m6anet("data.site_proba.csv",
#'                     probability_threshold = 0.9,
#'                     transcript_col = "transcript_id",
#'                     position_col   = "transcript_position",
#'                     prob_col       = "probability_modified")
#'
#' # Retain read depth for QC
#' m6a <- parse_m6anet("data.site_proba.csv", keep_extra = "n_reads")
#' }
#'
#' @import data.table
#'
#' @export
parse_m6anet <- function(m6anet_file,
                         probability_threshold = 0.9,
                         transcript_col = "transcript_id",
                         position_col   = "position",
                         prob_col       = "probability",
                         keep_extra     = NULL) {

  if (!file.exists(m6anet_file)) {
    stop("m6Anet file not found: ", m6anet_file)
  }
  if (!is.numeric(probability_threshold) ||
      probability_threshold < 0 || probability_threshold > 1) {
    stop("probability_threshold must be between 0 and 1")
  }

  m6a_data <- data.table::fread(m6anet_file)

  if (nrow(m6a_data) == 0) {
    stop("m6Anet file is empty: ", m6anet_file)
  }

  # Standardise column names
  data.table::setnames(
    m6a_data,
    old = c(transcript_col, position_col, prob_col),
    new = c("transcript_id", "position", "probability"),
    skip_absent = TRUE
  )

  required_cols <- c("transcript_id", "position", "probability")
  missing_cols  <- setdiff(required_cols, names(m6a_data))
  if (length(missing_cols) > 0) {
    stop("m6Anet output is missing columns: ", paste(missing_cols, collapse = ", "),
         "\n  Found: ", paste(names(m6a_data), collapse = ", "),
         "\n  Set transcript_col / position_col / prob_col to match your file.")
  }

  n_input <- nrow(m6a_data)

  # Drop NA probabilities
  na_prob_rows <- m6a_data[is.na(probability), .N]
  if (na_prob_rows > 0) {
    message("Dropping ", na_prob_rows, " row(s) with NA probability.")
    m6a_data <- m6a_data[!is.na(probability)]
  }

  # Filter by threshold
  m6a_data <- m6a_data[probability >= probability_threshold]

  if (nrow(m6a_data) == 0) {
    warning("No sites passed probability_threshold = ", probability_threshold,
            " (input had ", n_input, " rows).", call. = FALSE)
  }

  # m6Anet reports 0-based positions; convert to 1-based.
  # Verified against transcript sequence: only the shifted position gives A at
  # every site and reproduces the reported kmer.
  m6a_data[, position := as.integer(position) + 1L]

  # ── Column retention ────────────────────────────────────────────────────────
  # kmer is kept automatically for annotate_drach(). Anything named in
  # keep_extra is kept too - notably 'condition', which enables condition-level
  # comparison in annotate_m6a_switches_genomic().
  keep_cols <- c("transcript_id", "position", "probability")
  if ("kmer" %in% names(m6a_data)) {
    keep_cols <- c(keep_cols, "kmer")
  }
  if (!is.null(keep_extra)) {
    found   <- intersect(keep_extra, names(m6a_data))
    absent  <- setdiff(keep_extra, names(m6a_data))
    if (length(absent) > 0) {
      warning("keep_extra column(s) not found and ignored: ",
              paste(absent, collapse = ", "), call. = FALSE)
    }
    keep_cols <- unique(c(keep_cols, found))
  }
  m6a_data <- m6a_data[, keep_cols, with = FALSE]

  data.table::setorder(m6a_data, transcript_id, position)

  message(sprintf("Parsed %s: %d of %d sites retained (threshold %.2f)",
                  basename(m6anet_file), nrow(m6a_data), n_input,
                  probability_threshold))

  m6a_data
}
