#' Export Annotated m6A Switches
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param output_prefix Prefix for output files (default: "m6aswitch_results")
#' @param format Output format: "csv", "bed", "both" (default: "both")
#' @param add_igv_track Logical: generate IGV-compatible BED file (default: TRUE)
#'
#' @return List with paths to output files
#'
#' @examples
#' \dontrun{
#' export_annotated_switches(m6a_switches, output_prefix = "glioma_m6a")
#' }
#'
#' @import data.table
#'
#' @export
export_annotated_switches <- function(m6a_switches,
                                      output_prefix = "m6aswitch_results",
                                      format = "both",
                                      add_igv_track = TRUE) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  format <- match.arg(format, c("csv", "bed", "both"))
  output_files <- list()

  if (format %in% c("csv", "both")) {
    csv_file <- sprintf("%s.csv", output_prefix)
    data.table::fwrite(m6a_switches, file = csv_file)
    output_files$csv <- csv_file
    message("Wrote CSV: ", csv_file)
  }

  if (format %in% c("bed", "both") && add_igv_track) {
    bed_file <- sprintf("%s.bed", output_prefix)

    has_genomic <- "start" %in% names(m6a_switches)

    if (has_genomic) {
      chrom_col  <- m6a_switches$seqname
      start_col  <- m6a_switches$start - 1L
      end_col    <- m6a_switches$start
      strand_col <- m6a_switches$strand
    } else {
      chrom_col  <- m6a_switches$gene_id
      start_col  <- m6a_switches$position - 1L
      end_col    <- m6a_switches$position
      strand_col <- "+"
    }

    prob_score <- ifelse(
      is.na(m6a_switches$probability_a),
      m6a_switches$probability_b,
      m6a_switches$probability_a
    )

    bed_data <- data.table::data.table(
      chrom      = chrom_col,
      chromStart = start_col,
      chromEnd   = end_col,
      name       = m6a_switches$m6a_fate,
      score      = round(prob_score * 100, 0),
      strand     = strand_col,
      itemRGB    = ifelse(m6a_switches$m6a_fate == "LOST",   "255,0,0",
                   ifelse(m6a_switches$m6a_fate == "GAINED", "0,255,0", "0,0,255"))
    )

    data.table::setorder(bed_data, chrom, chromStart)

    con <- file(bed_file, "w")
    writeLines("# m6A switches - color: LOST=red, GAINED=green, RETAINED=blue", con)
    close(con)

    data.table::fwrite(bed_data, file = bed_file, sep = "\t",
                       col.names = TRUE, append = TRUE)
    output_files$bed <- bed_file
    message("Wrote BED: ", bed_file)
  }

  return(output_files)
}

#' Generate Summary Report
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param output_file Path to write report (default: "m6aswitch_report.txt")
#'
#' @return invisibly returns character vector of report lines
#'
#' @import data.table
#'
#' @export
generate_summary_report <- function(m6a_switches, output_file = "m6aswitch_report.txt") {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  report <- c(
    "=== m6Aswitch Analysis Report ===",
    sprintf("Analysis timestamp: %s", Sys.time()),
    "",
    "## Summary Statistics ##",
    sprintf("Total isoform switches analyzed: %d",
            data.table::uniqueN(m6a_switches[, .(gene_id, isoform_a, isoform_b)])),
    sprintf("Total m6A sites detected: %d", nrow(m6a_switches)),
    sprintf("  - LOST: %d (%.1f%%)",
            nrow(m6a_switches[m6a_fate == "LOST"]),
            nrow(m6a_switches[m6a_fate == "LOST"]) / nrow(m6a_switches) * 100),
    sprintf("  - GAINED: %d (%.1f%%)",
            nrow(m6a_switches[m6a_fate == "GAINED"]),
            nrow(m6a_switches[m6a_fate == "GAINED"]) / nrow(m6a_switches) * 100),
    sprintf("  - RETAINED: %d (%.1f%%)",
            nrow(m6a_switches[m6a_fate == "RETAINED"]),
            nrow(m6a_switches[m6a_fate == "RETAINED"]) / nrow(m6a_switches) * 100),
    "",
    "## Top Genes with m6A Changes ##",
    ""
  )

  top_genes <- m6a_switches[, .N, by = gene_id][order(-N)][1:min(10, nrow(m6a_switches))]
  for (i in seq_len(nrow(top_genes))) {
    report <- c(report, sprintf("  %s: %d m6A events",
                                top_genes[i, gene_id], top_genes[i, N]))
  }

  report <- c(report, "")
  writeLines(report, con = output_file)
  message("Wrote report: ", output_file)

  invisible(report)
}
