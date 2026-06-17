#' Annotate m6A Changes Across Isoform Switches
#'
#' Core function: integrates m6A sites and isoform switches to identify which 
#' m6A sites are gained or lost during isoform switching events.
#'
#' @param m6a_sites data.table from parse_m6anet()
#' @param iso_switches data.table from parse_isoform_switch()
#' @param iso_sequences data.table with columns: isoform_id, sequence
#'        (transcript sequences; required to map m6A positions)
#'
#' @return data.table with columns:
#'   - gene_id, isoform_a, isoform_b (switch info)
#'   - position (m6A position in genomic coords)
#'   - m6a_in_isoform_a, m6a_in_isoform_b (TRUE/FALSE)
#'   - m6a_fate (LOST, GAINED, RETAINED)
#'   - probability_a, probability_b (m6A probabilities in each isoform)
#'   - fdr (isoform switch FDR)
#'
#' @examples
#' \dontrun{
#' m6a_switches <- annotate_m6a_switches(m6a_sites, iso_switches, iso_sequences)
#' }
#'
#' @import data.table
#' @import GenomicRanges
#'
#' @export
annotate_m6a_switches <- function(m6a_sites, iso_switches, iso_sequences) {
  
  # Validate inputs
  if (!is.data.table(m6a_sites)) {
    stop("m6a_sites must be a data.table from parse_m6anet()")
  }
  if (!is.data.table(iso_switches)) {
    stop("iso_switches must be a data.table from parse_isoform_switch()")
  }
  if (!is.data.table(iso_sequences)) {
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
  
  for (i in 1:nrow(iso_switches)) {
    switch_row <- iso_switches[i]
    gene <- switch_row$gene_id
    iso_a <- switch_row$isoform_a
    iso_b <- switch_row$isoform_b
    fdr <- switch_row$fdr
    
    # Get m6A sites for both isoforms
    m6a_in_a <- iso_m6a_map[isoform_id == iso_a]
    m6a_in_b <- iso_m6a_map[isoform_id == iso_b]
    
    # Create combinations
    if (nrow(m6a_in_a) > 0 | nrow(m6a_in_b) > 0) {
      all_positions <- unique(c(m6a_in_a$position, m6a_in_b$position))
      
      for (pos in all_positions) {
        prob_a <- m6a_in_a[position == pos, probability]
        prob_b <- m6a_in_b[position == pos, probability]
        
        in_a <- length(prob_a) > 0
        in_b <- length(prob_b) > 0
        
        # Determine fate
        if (in_a & !in_b) {
          fate <- "LOST"
        } else if (!in_a & in_b) {
          fate <- "GAINED"
        } else {
          fate <- "RETAINED"
        }
        
        result_list[[length(result_list) + 1]] <- data.table(
          gene_id = gene,
          isoform_a = iso_a,
          isoform_b = iso_b,
          position = pos,
          m6a_in_isoform_a = in_a,
          m6a_in_isoform_b = in_b,
          m6a_fate = fate,
          probability_a = if (length(prob_a) > 0) prob_a[1] else NA_real_,
          probability_b = if (length(prob_b) > 0) prob_b[1] else NA_real_,
          fdr = fdr
        )
      }
    }
  }
  
  if (length(result_list) == 0) {
    warning("No m6A sites found in isoform switches. Returning empty result.")
    return(data.table())
  }
  
  result <- rbindlist(result_list)
  setorder(result, gene_id, fdr, m6a_fate)
  
  return(result)
}
