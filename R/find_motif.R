#' Find DRACH Motif at m6A Site
#'
#' Checks whether the base at `position` is the central A of a DRACH motif.
#' DRACH = D-R-A-C-H where:
#'   D = A/G/U, R = A/G, A = A, C = C, H = A/C/U
#'
#' @param sequence Character string of RNA sequence
#' @param position Numeric position within sequence (1-indexed)
#' @param context_bp Bases before/after position to check (kept for API compatibility;
#'   DRACH check uses a strict 5-mer centered on `position`)
#'
#' @return logical TRUE if DRACH motif centered at `position`, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' find_motif("AGACA", 3)  # TRUE
#' }
#'
#' @import Biostrings
#'
#' @export
find_motif <- function(sequence, position, context_bp = 2) {

  # Validate inputs
  if (!is.character(sequence) || length(sequence) != 1) {
    stop("sequence must be a single character string")
  }
  if (!is.numeric(position) || length(position) != 1 || is.na(position)) {
    stop("position must be a single non-NA numeric value")
  }

  position <- as.integer(position)
  if (position < 1 || position > nchar(sequence)) {
    stop("position must be numeric and within sequence bounds")
  }

  # Convert DNA->RNA and uppercase
  sequence <- toupper(chartr("T", "U", sequence))

  # Need full 5-mer around center
  if (position - 2 < 1 || position + 2 > nchar(sequence)) {
    return(FALSE)
  }

  motif <- substr(sequence, position - 2, position + 2)
  # DRACH centered at A: [AGU][AG]AC[ACU]
  return(grepl("^[AGU][AG]AC[ACU]$", motif))
}

#' Annotate m6A Sites with DRACH Motif Status
#'
#' Adds DRACH motif annotation to m6A switching results.
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param sequences data.table with columns: isoform_id, sequence
#'
#' @return m6a_switches with added columns:
#'   - drach_motif_a (logical)
#'   - drach_motif_b (logical)
#'   - drach_motif (logical summary)
#'
#' @import data.table
#'
#' @export
annotate_drach <- function(m6a_switches, sequences) {

  if (!data.table::is.data.table(m6a_switches)) {
    stop("m6a_switches must be a data.table")
  }
  if (!data.table::is.data.table(sequences)) {
    stop("sequences must be a data.table")
  }
  if (!all(c("isoform_id", "sequence") %in% names(sequences))) {
    stop("sequences must contain columns: isoform_id, sequence")
  }

  # Handle empty input
  if (nrow(m6a_switches) == 0) {
    m6a_switches[, drach_motif_a := logical(0)]
    m6a_switches[, drach_motif_b := logical(0)]
    m6a_switches[, drach_motif := logical(0)]
    return(m6a_switches)
  }

  # Reset motif columns each run (prevents old type contamination)
  m6a_switches[, c("drach_motif_a", "drach_motif_b", "drach_motif") := NULL]
  m6a_switches[, drach_motif_a := as.logical(NA)]
  m6a_switches[, drach_motif_b := as.logical(NA)]

  seq_map <- unique(sequences[, .(isoform_id, sequence)], by = "isoform_id")

  for (i in seq_len(nrow(m6a_switches))) {
    iso_a <- m6a_switches[i, isoform_a]
    iso_b <- m6a_switches[i, isoform_b]
    pos  <- m6a_switches[i, position]

    in_a <- isTRUE(m6a_switches[i, m6a_in_isoform_a])
    in_b <- isTRUE(m6a_switches[i, m6a_in_isoform_b])

    if (is.na(pos)) next

    if (in_a) {
      seq_a <- seq_map[isoform_id == iso_a, sequence]
      if (length(seq_a) > 0 && !is.na(seq_a[1])) {
        m6a_switches[i, drach_motif_a := tryCatch(
          find_motif(seq_a[1], pos),
          error = function(e) as.logical(NA)
        )]
      }
    }

    if (in_b) {
      seq_b <- seq_map[isoform_id == iso_b, sequence]
      if (length(seq_b) > 0 && !is.na(seq_b[1])) {
        m6a_switches[i, drach_motif_b := tryCatch(
          find_motif(seq_b[1], pos),
          error = function(e) as.logical(NA)
        )]
      }
    }
  }

  # Force logical type before summary
  m6a_switches[, drach_motif_a := as.logical(drach_motif_a)]
  m6a_switches[, drach_motif_b := as.logical(drach_motif_b)]

  # Row-wise summary
  m6a_switches[, drach_motif := data.table::fifelse(
    (isTRUE(drach_motif_a) || isTRUE(drach_motif_b)), TRUE,
    data.table::fifelse(
      (!is.na(drach_motif_a) && !is.na(drach_motif_b) &&
         identical(drach_motif_a, FALSE) && identical(drach_motif_b, FALSE)),
      FALSE,
      as.logical(NA)
    )
  ), by = seq_len(nrow(m6a_switches))]

  return(m6a_switches)
}
