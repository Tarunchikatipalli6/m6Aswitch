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
  
  # DRACH is typically centered on A
  # Search within window for pattern: [AGAU]RACH or [AGAU]ACH
  
  drach_pattern <- "[AGAU]R[AU]CH"  # Simplified: D=A/G/U, R=A/G, A=A/U, C=C, H=A/C/U
  
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
  
  # Merge sequences
  seq_map <- sequences[, .(isoform_id, sequence)]
  
  m6a_switches[, drach_motif := NA]
  
  for (i in 1:nrow(m6a_switches)) {
    iso_a <- m6a_switches[i, isoform_a]
    iso_b <- m6a_switches[i, isoform_b]
    position <- m6a_switches[i, position]
    
    # Check DRACH in isoform A (if m6A present)
    if (m6a_switches[i, m6a_in_isoform_a]) {
      seq_a <- seq_map[isoform_id == iso_a, sequence]
      if (length(seq_a) > 0) {
        m6a_switches[i, drach_motif_a := find_motif(seq_a, position)]
      }
    }
    
    # Check DRACH in isoform B (if m6A present)
    if (m6a_switches[i, m6a_in_isoform_b]) {
      seq_b <- seq_map[isoform_id == iso_b, sequence]
      if (length(seq_b) > 0) {
        m6a_switches[i, drach_motif_b := find_motif(seq_b, position)]
      }
    }
  }
  
  return(m6a_switches)
}
