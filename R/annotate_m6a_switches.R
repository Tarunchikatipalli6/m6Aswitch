#' Annotate m6A Sites with DRACH Motif Status
#'
#' Adds DRACH motif annotation to m6A switching results using the kmer
#' sequence context already provided by m6Anet output. This is simpler
#' and more accurate than computing the motif from transcript sequences.
#'
#' DRACH = D-R-A-C-H where:
#'   D = A/G/U, R = A/G, A = A, C = C, H = A/C/U
#'
#' @param m6a_switches data.table from annotate_m6a_switches_genomic()
#' @param m6a_condition_a data.table from parse_m6anet() for condition A
#'   (must contain kmer column — requires m6Anet output with kmer)
#' @param m6a_condition_b data.table from parse_m6anet() for condition B
#'   (must contain kmer column — requires m6Anet output with kmer)
#'
#' @return m6a_switches with added column:
#'   - drach_motif (logical): TRUE if the m6A site is in a DRACH motif
#'
#' @details
#' m6Anet outputs a kmer column containing the 5-mer sequence context
#' centered on each m6A site. This function uses that kmer directly
#' to check the DRACH pattern, removing the need for transcript
#' sequence files entirely.
#'
#' If kmer column is not present in m6A data, the function returns
#' drach_motif as NA for all rows with a warning.
#'
#' @examples
#' \dontrun{
#' m6a_final <- annotate_drach(m6a_switches, m6a_igf, m6a_naive)
#' table(m6a_final$drach_motif)
#' }
#'
#' @import data.table
#'
#' @export
annotate_drach <- function(m6a_switches, m6a_condition_a, m6a_condition_b) {

  if (!data.table::is.data.table(m6a_switches)) {
    stop("m6a_switches must be a data.table from annotate_m6a_switches_genomic()")
  }
  if (!data.table::is.data.table(m6a_condition_a)) {
    stop("m6a_condition_a must be a data.table from parse_m6anet()")
  }
  if (!data.table::is.data.table(m6a_condition_b)) {
    stop("m6a_condition_b must be a data.table from parse_m6anet()")
  }

  # Handle empty input
  if (nrow(m6a_switches) == 0) {
    m6a_switches[, drach_motif := logical(0)]
    return(m6a_switches)
  }

  # Check kmer column exists
  if (!"kmer" %in% names(m6a_condition_a) ||
      !"kmer" %in% names(m6a_condition_b)) {
    warning("kmer column not found in m6A data. drach_motif will be NA. ",
            "Make sure parse_m6anet() is up to date.")
    m6a_switches[, drach_motif := NA]
    return(m6a_switches)
  }

  # Build kmer lookup from both conditions combined
  kmer_map <- unique(rbind(
    m6a_condition_a[, .(transcript_id, position, kmer)],
    m6a_condition_b[, .(transcript_id, position, kmer)]
  ), by = c("transcript_id", "position"))

  # Determine position column
  pos_col <- if ("transcript_position" %in% names(m6a_switches)) {
    "transcript_position"
  } else {
    "position"
  }

  # Join kmer to switches using isoform_a transcript and position
  # isoform_a is the reference isoform for position lookup
  m6a_switches <- merge(
    m6a_switches,
    kmer_map[, .(transcript_id, position, kmer)],
    by.x = c("isoform_a", pos_col),
    by.y = c("transcript_id", "position"),
    all.x = TRUE
  )

  # Check DRACH pattern on kmer
  # Convert T to U (DNA to RNA) and check [AGU][AG]AC[ACU]
  m6a_switches[, drach_motif := data.table::fifelse(
    !is.na(kmer),
    grepl(
      "^[AGU][AG]AC[ACU]$",
      toupper(chartr("T", "U", kmer))
    ),
    NA
  )]

  # Remove kmer column from output — internal use only
  m6a_switches[, kmer := NULL]

  return(m6a_switches)
}
