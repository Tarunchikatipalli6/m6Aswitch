#' Find DRACH Motif at m6A Site
#'
#' Checks if an m6A site location matches a DRACH motif pattern 
#' (D=A/G/U, R=A/G) in the surrounding sequence context.
#'
#' @param sequence Character string of RNA sequence
#' @param position Numeric position within sequence (1-indexed)
#' @param context_bp Bases before/after position to check (default: 2)
#'
#' @return logical TRUE if DRACH motif present, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' find_motif("ACUGGACC", 4)  # Check if position 4 (A) is in DRACH
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
  if (!is.numeric(position) || position < 1 || position > nchar(sequence)) {
    stop("position must be numeric and within sequence bounds")
  }
  
  # Convert to uppercase
  sequence <- toupper(sequence)
  
  # Extract context window
  start <- max(1, position - context_bp)
  end <- min(nchar(sequence), position + context_bp)
  window <- substr(sequence, start, end)
  
  # Check for DRACH motif pattern
  # D = A/G/U (any purine or U)
  # R = A/G (purine)
  # A = adenosine (target for m6A)
  # C = cytosine
  # H = A/C/U (any except G)
  
  # DRACH pattern: [AGU][AG][AU]C[AUC]
  # More flexible: look for D-R-A-C-H where:
  # - D is at any position in window
  # - R follows D
  # - A follows R
  # - C follows A
  # - H follows C
  
  drach_pattern <- "[AGU][AG][AU]C[AUC]"
  
  if (grepl(drach_pattern, window, ignore.case = TRUE)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

#' Annotate m6A Sites with DRACH Motif Status
#'
#' Adds DRACH motif annotation to m6A switching results.
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param sequences data.table with columns: isoform_id, sequence
#'
#' @return m6a_switches with added column: drach_motif (logical)
#'
#' @import data.table
#'
#' @export
annotate_drach <- function(m6a_switches, sequences) {
  
  if (!is.data.table(m6a_switches)) {
    stop("m6a_switches must be a data.table")
  }
  
  # Handle empty input
  if (nrow(m6a_switches) == 0) {
    m6a_switches[, drach_motif := logical(0)]
    return(m6a_switches)
  }
  
  # Merge sequences
  seq_map <- sequences[, .(isoform_id, sequence)]
  
  m6a_switches[, drach_motif := NA]
  m6a_switches[, drach_motif_a := NA_character_]
  m6a_switches[, drach_motif_b := NA_character_]
  
  for (i in 1:nrow(m6a_switches)) {
    iso_a <- m6a_switches[i, isoform_a]
    iso_b <- m6a_switches[i, isoform_b]
    position <- m6a_switches[i, position]
    m6a_in_a <- m6a_switches[i, m6a_in_isoform_a]
    m6a_in_b <- m6a_switches[i, m6a_in_isoform_b]
    
    # Skip if position is NA or not numeric
    if (is.na(position) || !is.numeric(position)) {
      next
    }
    
    # Check DRACH in isoform A (if m6A present)
    if (!is.na(m6a_in_a) && m6a_in_a) {
      seq_a_vec <- seq_map[isoform_id == iso_a, sequence]
      if (length(seq_a_vec) > 0) {
        seq_a <- seq_a_vec[1]  # Extract first element
        tryCatch({
          m6a_switches[i, drach_motif_a := find_motif(seq_a, position)]
        }, error = function(e) {
          NULL
        })
      }
    }
    
    # Check DRACH in isoform B (if m6A present)
    if (!is.na(m6a_in_b) && m6a_in_b) {
      seq_b_vec <- seq_map[isoform_id == iso_b, sequence]
      if (length(seq_b_vec) > 0) {
        seq_b <- seq_b_vec[1]  # Extract first element
        tryCatch({
          m6a_switches[i, drach_motif_b := find_motif(seq_b, position)]
        }, error = function(e) {
          NULL
        })
      }
    }
  }
  
  return(m6a_switches)
}
