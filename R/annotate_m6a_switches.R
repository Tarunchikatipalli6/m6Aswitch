#' Annotate m6A Changes Across Isoform Switches (Transcript Coordinates)
#'
#' DEPRECATED FOR PUBLICATION USE.
#' 
#' Core function: integrates m6A sites and isoform switches to identify which
#' m6A sites are gained or lost during isoform switching events. 
#' 
#' WARNING: This function compares transcript-level positions directly, which can
#' produce FALSE calls when isoforms have different exon structures (see Details).
#' For publication-quality results, use annotate_m6a_switches_genomic() instead.
#'
#' @param m6a_sites data.table from parse_m6anet()
#' @param iso_switches data.table from parse_isoform_switch()
#' @param iso_sequences data.table with columns: isoform_id, sequence
#'
#' @return data.table with columns:
#'   - gene_id, isoform_a, isoform_b (switch info)
#'   - condition_1, condition_2 (condition labels from the comparison)
#'   - position (m6A position in TRANSCRIPT coordinates - see Details)
#'   - m6a_in_isoform_a, m6a_in_isoform_b (TRUE/FALSE)
#'   - m6a_fate (LOST, GAINED, RETAINED)
#'   - probability_a, probability_b (m6A probabilities in each isoform)
#'   - fdr, dif (isoform switch FDR and direction indicator)
#'
#' @details
#' CRITICAL CAVEAT: This function operates on transcript-level coordinates.
#' Position 150 in isoform A and position 150 in isoform B are counted from the
#' start of each transcript independently. If the isoforms have different exon
#' structures (e.g., one skips an exon), these identical position numbers may
#' correspond to completely different genomic locations.
#'
#' Example of the problem:
#'   Isoform A: Exon1 (chr10:1000-1100) + Exon2 (chr10:5000-5200)
#'   Isoform B: Exon1 (chr10:1000-1100) + Exon3 (chr10:8000-8200)  [Exon2 skipped]
#'
#'   Position 150 in both transcripts:
#'   - Isoform A position 150 = genomic chr10:5050 (in Exon2)
#'   - Isoform B position 150 = genomic chr10:8050 (in Exon3)
#'
#'   This function would compare these as "same position" and report RETAINED,
#'   but they are actually at different genomic loci.
#'
#' FOR PUBLICATION: Use annotate_m6a_switches_genomic() with a GTF file instead.
#' This lifts coordinates to the genome and enables scientifically valid comparisons.
#'
#' m6A fate categories (transcript-level, use with caution):
#' - LOST: m6A at this transcript position in isoform_a but not isoform_b
#' - GAINED: m6A at this transcript position in isoform_b but not isoform_a
#' - RETAINED: m6A at this transcript position in both isoforms
#'
#' @examples
#' \dontrun{
#' # NOT RECOMMENDED FOR PUBLICATION
#' # Use this only for exploratory analysis or when isoforms share identical exon structure
#' m6a_switches <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
#' }
#'
#' @seealso annotate_m6a_switches_genomic() for publication-quality genomic coordinate analysis
#'
#' @import data.table
#'
#' @export
annotate_m6a_switches <- function(m6a_sites, iso_switches, iso_sequences) {

  # Validate inputs
  if (!data.table::is.data.table(m6a_sites)) {
    stop("m6a_sites must be a data.table from parse_m6anet()")
  }
  if (!data.table::is.data.table(iso_switches)) {
    stop("iso_switches must be a data.table from parse_isoform_switch()")
  }
  if (!data.table::is.data.table(iso_sequences)) {
    stop("iso_sequences must be a data.table with columns: isoform_id, sequence")
  }

  # Step 1: Build mapping of isoform -> m6A sites
  iso_m6a_map <- m6a_sites[, .(
    isoform_id = transcript_id,
    position = position,
    probability = probability
  )]

  # Step 2: For each isoform switch, check which m6A sites are present in A vs B
  result_list <- list()

  if (nrow(iso_switches) == 0) {
    warning("iso_switches is empty. Returning empty result.")
    return(data.table())
  }

  for (i in seq_len(nrow(iso_switches))) {
    switch_row <- iso_switches[i]
    gene <- switch_row$gene_id
    iso_a <- switch_row$isoform_a
    iso_b <- switch_row$isoform_b
    fdr <- switch_row$fdr

    # Extract condition information (pairwise comparison)
    condition_1 <- if ("condition_1" %in% names(switch_row))
      switch_row$condition_1 else NA_character_
    condition_2 <- if ("condition_2" %in% names(switch_row))
      switch_row$condition_2 else NA_character_

    # Extract dIF (direction indicator)
    dif <- if ("dif" %in% names(switch_row))
      switch_row$dif else NA_real_

    # Get m6A sites for both isoforms
    m6a_in_a <- iso_m6a_map[isoform_id == iso_a]
    m6a_in_b <- iso_m6a_map[isoform_id == iso_b]

    # Create combinations
    if (nrow(m6a_in_a) > 0 || nrow(m6a_in_b) > 0) {
      all_positions <- unique(c(m6a_in_a$position, m6a_in_b$position))

      for (pos in all_positions) {
        prob_a <- m6a_in_a[position == pos, probability]
        prob_b <- m6a_in_b[position == pos, probability]

        in_a <- length(prob_a) > 0
        in_b <- length(prob_b) > 0

        # Determine fate
        if (in_a && !in_b) {
          fate <- "LOST"
        } else if (!in_a && in_b) {
          fate <- "GAINED"
        } else {
          fate <- "RETAINED"
        }

        result_list[[length(result_list) + 1]] <- data.table(
          gene_id = gene,
          isoform_a = iso_a,
          isoform_b = iso_b,
          condition_1 = condition_1,
          condition_2 = condition_2,
          position = pos,
          m6a_in_isoform_a = in_a,
          m6a_in_isoform_b = in_b,
          m6a_fate = fate,
          probability_a = if (length(prob_a) > 0) prob_a[1] else NA_real_,
          probability_b = if (length(prob_b) > 0) prob_b[1] else NA_real_,
          fdr = fdr,
          dif = dif
        )
      }
    }
  }

  if (length(result_list) == 0) {
    warning("No m6A sites found in isoform switches. Returning empty result.")
    return(data.table())
  }

  result <- data.table::rbindlist(result_list)
  data.table::setorder(result, gene_id, fdr, m6a_fate)

  return(result)
}

#' Find DRACH Motif at m6A Site
#'
#' Checks whether the base at `position` is the central A of a DRACH motif.
#' DRACH = D-R-A-C-H where:
#'   D = A/G/U, R = A/G, A = A, C = C, H = A/C/U
#'
#' @param sequence Character string of RNA sequence
#' @param position Numeric position within sequence (1-indexed)
#' @param context_bp Bases before/after position to check (default: 2)
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
    stop("position must be within sequence bounds")
  }

  # Convert DNA->RNA and uppercase
  sequence <- toupper(chartr("T", "U", sequence))

  # Need full 5-mer centered at position (2 upstream + center + 2 downstream)
  if (position - 2 < 1 || position + 2 > nchar(sequence)) {
    return(FALSE)
  }

  motif <- substr(sequence, position - 2, position + 2)
  # DRACH centered at A: [AGU][AG]AC[ACU]
  return(grepl("^[AGU][AG]AC[ACU]$", motif))
}

#' Annotate m6A Sites with DRACH Motif Status
#'
#' Adds DRACH motif annotations to m6A switching results.
#'
#' @param m6a_switches data.table from annotate_m6a_switches() or 
#'        annotate_m6a_switches_genomic()
#' @param sequences data.table with columns: isoform_id, sequence
#'
#' @return m6a_switches with added columns:
#'   - drach_motif_a (logical)
#'   - drach_motif_b (logical)
#'   - drach_motif (logical summary; TRUE if either side is TRUE,
#'     FALSE if both present and FALSE, otherwise NA)
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

  # Keep first sequence per isoform_id
  seq_map <- unique(sequences[, .(isoform_id, sequence)], by = "isoform_id")

  m6a_switches[, drach_motif_a := as.logical(NA)]
  m6a_switches[, drach_motif_b := as.logical(NA)]

  # Determine which column contains position information
  if ("position" %in% names(m6a_switches)) {
    pos_col <- "position"
  } else if ("start" %in% names(m6a_switches)) {
    pos_col <- "start"  # Genomic coordinates
  } else {
    stop("m6a_switches must contain either 'position' or 'start' column")
  }

  for (i in seq_len(nrow(m6a_switches))) {
    iso_a <- m6a_switches[i, isoform_a]
    iso_b <- m6a_switches[i, isoform_b]
    pos <- m6a_switches[i, get(pos_col)]
    in_a <- isTRUE(m6a_switches[i, m6a_in_isoform_a])
    in_b <- isTRUE(m6a_switches[i, m6a_in_isoform_b])

    if (is.na(pos)) next

    # Isoform A
    if (in_a) {
      seq_a <- seq_map[isoform_id == iso_a, sequence]
      if (length(seq_a) > 0 && !is.na(seq_a[1])) {
        m6a_switches[i, drach_motif_a := tryCatch(
          find_motif(seq_a[1], pos),
          error = function(e) as.logical(NA)
        )]
      }
    }

    # Isoform B
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

  # Summary column
  m6a_switches[, drach_motif := fifelse(
    isTRUE(drach_motif_a) | isTRUE(drach_motif_b), TRUE,
    fifelse(!is.na(drach_motif_a) & !is.na(drach_motif_b) &
              !drach_motif_a & !drach_motif_b, FALSE, as.logical(NA))
  ), by = seq_len(nrow(m6a_switches))]

  return(m6a_switches)
}
