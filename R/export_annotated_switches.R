#' Export Annotated m6A Switches
#'
#' Writes results to CSV and/or BED. The BED \code{name} and colour fields use
#' whichever classification you select via \code{color_by}.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#' @param output_prefix Prefix for output files (default "m6aswitch_results").
#' @param format One of "csv", "bed", "both" (default "both").
#' @param add_igv_track Logical: write an IGV track line at the top of the BED
#'   (default TRUE).
#' @param color_by Which classification to encode in the BED: \code{"isoform_status"}
#'   (structural, default) or \code{"m6a_fate"} (regulatory; requires condition
#'   information supplied upstream).
#'
#' @return List with paths to output files.
#'
#' @section BED colours:
#' \code{isoform_status}: ISOFORM_A_ONLY red, ISOFORM_B_ONLY green,
#' IN_BOTH_ISOFORMS blue.
#'
#' \code{m6a_fate}: LOST red, GAINED green, RETAINED blue.
#'
#' @examples
#' \dontrun{
#' export_annotated_switches(res, output_prefix = "glioma_m6a")
#' export_annotated_switches(res, "glioma_cond", color_by = "m6a_fate")
#' }
#'
#' @import data.table
#'
#' @export
export_annotated_switches <- function(m6a_switches,
                                      output_prefix = "m6aswitch_results",
                                      format        = "both",
                                      add_igv_track = TRUE,
                                      color_by      = c("isoform_status", "m6a_fate")) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  format   <- match.arg(format, c("csv", "bed", "both"))
  color_by <- match.arg(color_by)

  output_files <- list()

  if (format %in% c("csv", "both")) {
    csv_file <- sprintf("%s.csv", output_prefix)
    data.table::fwrite(m6a_switches, file = csv_file)
    output_files$csv <- csv_file
    message("Wrote CSV: ", csv_file)
  }

  if (format %in% c("bed", "both")) {

    if (!color_by %in% names(m6a_switches)) {
      stop("Column '", color_by, "' not found. Cannot write BED. Available: ",
           paste(names(m6a_switches), collapse = ", "))
    }
    if (all(is.na(m6a_switches[[color_by]]))) {
      stop("Column '", color_by, "' is entirely NA. ",
           if (color_by == "m6a_fate")
             "Condition information was not supplied upstream; use color_by = 'isoform_status'."
           else "")
    }

    bed_file <- sprintf("%s.bed", output_prefix)

    has_genomic <- all(c("seqname", "start", "strand") %in% names(m6a_switches))

    if (has_genomic) {
      chrom_col  <- m6a_switches$seqname
      start_col  <- m6a_switches$start - 1L    # BED is 0-based half-open
      end_col    <- if ("end" %in% names(m6a_switches)) {
        m6a_switches$end
      } else {
        m6a_switches$start
      }
      strand_col <- m6a_switches$strand
    } else {
      warning("No genomic coordinates found; using gene_id as the BED chrom ",
              "field. The result will not load in a genome browser.",
              call. = FALSE)
      chrom_col  <- m6a_switches$gene_id
      start_col  <- m6a_switches$position - 1L
      end_col    <- m6a_switches$position
      strand_col <- "+"
    }

    # Score: probability from whichever isoform carries the site
    prob_score <- data.table::fcoalesce(
      m6a_switches$probability_a,
      m6a_switches$probability_b
    )
    score_col <- ifelse(is.na(prob_score), 0L,
                        pmin(1000L, round(prob_score * 1000)))

    # Colour lookup for the selected classification
    rgb_map <- if (color_by == "isoform_status") {
      c(ISOFORM_A_ONLY   = "255,0,0",
        ISOFORM_B_ONLY   = "0,255,0",
        IN_BOTH_ISOFORMS = "0,0,255")
    } else {
      c(LOST     = "255,0,0",
        GAINED   = "0,255,0",
        RETAINED = "0,0,255")
    }

    class_vals <- as.character(m6a_switches[[color_by]])
    item_rgb   <- unname(rgb_map[class_vals])
    item_rgb[is.na(item_rgb)] <- "128,128,128"    # unknown or NA

    bed_data <- data.table::data.table(
      chrom      = chrom_col,
      chromStart = as.integer(start_col),
      chromEnd   = as.integer(end_col),
      name       = ifelse(is.na(class_vals), "UNCLASSIFIED", class_vals),
      score      = as.integer(score_col),
      strand     = strand_col,
      thickStart = as.integer(start_col),
      thickEnd   = as.integer(end_col),
      itemRgb    = item_rgb
    )

    data.table::setorder(bed_data, chrom, chromStart)

    # Write the track line first, then append the records.
    # fwrite(append = TRUE) requires the file to exist, so create it here even
    # when no track line is requested.
    if (isTRUE(add_igv_track)) {
      writeLines(
        sprintf('track name="m6A %s" itemRgb="On"',
                if (color_by == "isoform_status") "isoform status" else "condition fate"),
        con = bed_file
      )
    } else {
      file.create(bed_file)
    }

    data.table::fwrite(bed_data, file = bed_file, sep = "\t",
                       col.names = FALSE, append = TRUE)

    output_files$bed <- bed_file
    message("Wrote BED: ", bed_file, "  (coloured by ", color_by, ")")
  }

  output_files
}


#' Generate Summary Report
#'
#' Writes a plain-text summary of an m6Aswitch analysis.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()},
#'   optionally after \code{classify_lost_mechanism()}.
#' @param output_file Path to write (default "m6aswitch_report.txt").
#'
#' @return Invisibly, a character vector of report lines.
#'
#' @import data.table
#' @importFrom utils head
#'
#' @export
generate_summary_report <- function(m6a_switches,
                                    output_file = "m6aswitch_report.txt") {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  n_total <- nrow(m6a_switches)

  report <- c(
    "=== m6Aswitch Analysis Report ===",
    sprintf("Generated: %s", Sys.time()),
    "",
    "## Overview ##",
    sprintf("Isoform switch pairs analysed : %d",
            data.table::uniqueN(m6a_switches[, .(gene_id, isoform_a, isoform_b)])),
    sprintf("Genes                         : %d",
            data.table::uniqueN(m6a_switches$gene_id)),
    sprintf("m6A sites                     : %d", n_total),
    ""
  )

  # ── Structural: which isoform carries the site ──────────────────────────────
  if ("isoform_status" %in% names(m6a_switches)) {
    report <- c(report,
      "## Isoform status (structural) ##",
      "Which transcript of the switch pair carries each site.",
      "")
    tab <- m6a_switches[, .N, by = isoform_status][order(-N)]
    for (i in seq_len(nrow(tab))) {
      report <- c(report, sprintf("  %-20s : %6d (%5.1f%%)",
                                  tab$isoform_status[i], tab$N[i],
                                  100 * tab$N[i] / n_total))
    }
    report <- c(report, "")
  }

  # ── Regulatory: which condition detected the site ───────────────────────────
  if ("m6a_fate" %in% names(m6a_switches)) {
    if (all(is.na(m6a_switches$m6a_fate))) {
      report <- c(report,
        "## Condition fate (regulatory) ##",
        "  Not computed - no condition information was supplied to",
        "  lift_m6a_to_genomic(). See ?annotate_m6a_switches_genomic.",
        "")
    } else {
      cond1 <- if ("condition_1" %in% names(m6a_switches)) m6a_switches$condition_1[1] else "condition 1"
      cond2 <- if ("condition_2" %in% names(m6a_switches)) m6a_switches$condition_2[1] else "condition 2"
      report <- c(report,
        "## Condition fate (regulatory) ##",
        sprintf("  LOST     = detected in %s only", cond1),
        sprintf("  GAINED   = detected in %s only", cond2),
        "  RETAINED = detected in both",
        "")
      tab <- m6a_switches[!is.na(m6a_fate), .N, by = m6a_fate][order(-N)]
      n_c <- sum(tab$N)
      for (i in seq_len(nrow(tab))) {
        report <- c(report, sprintf("  %-20s : %6d (%5.1f%%)",
                                    tab$m6a_fate[i], tab$N[i],
                                    100 * tab$N[i] / n_c))
      }
      report <- c(report, "")
    }
  }

  # ── Mechanism ───────────────────────────────────────────────────────────────
  if ("isoform_mechanism" %in% names(m6a_switches) &&
      any(!is.na(m6a_switches$isoform_mechanism))) {
    report <- c(report,
      "## Mechanism ##",
      "Structural = the sequence is absent from one isoform.",
      "Regulatory = the sequence is present in both, modification on one only.",
      "")
    tab <- m6a_switches[!is.na(isoform_mechanism), .N,
                        by = isoform_mechanism][order(-N)]
    n_m <- sum(tab$N)
    for (i in seq_len(nrow(tab))) {
      report <- c(report, sprintf("  %-24s : %6d (%5.1f%%)",
                                  tab$isoform_mechanism[i], tab$N[i],
                                  100 * tab$N[i] / n_m))
    }

    n_struct <- sum(m6a_switches$isoform_mechanism %in%
                    c("A_ONLY_EXON_SKIPPED", "B_ONLY_EXON_INCLUDED"), na.rm = TRUE)
    n_regul  <- sum(m6a_switches$isoform_mechanism %in%
                    c("A_ONLY_UNMETHYLATED", "B_ONLY_NEW_METHYLATION"), na.rm = TRUE)
    if (n_struct + n_regul > 0) {
      report <- c(report, "",
        sprintf("  Explained by splicing     : %d (%.1f%%)",
                n_struct, 100 * n_struct / (n_struct + n_regul)),
        sprintf("  Not explained by splicing : %d (%.1f%%)",
                n_regul,  100 * n_regul  / (n_struct + n_regul)))
    }
    report <- c(report, "")
  }

  # ── DRACH ───────────────────────────────────────────────────────────────────
  if ("drach_motif" %in% names(m6a_switches)) {
    n_drach <- sum(m6a_switches$drach_motif, na.rm = TRUE)
    n_known <- sum(!is.na(m6a_switches$drach_motif))
    if (n_known > 0) {
      report <- c(report,
        "## DRACH motif ##",
        sprintf("  In DRACH context: %d of %d (%.1f%%)",
                n_drach, n_known, 100 * n_drach / n_known),
        "  Note: m6Anet only evaluates DRACH positions, so 100% is expected",
        "  for m6Anet input.",
        "")
    }
  }

  # ── Top genes ───────────────────────────────────────────────────────────────
  report <- c(report, "## Top genes by m6A site count ##", "")
  top_genes <- utils::head(m6a_switches[, .N, by = gene_id][order(-N)], 10L)
  for (i in seq_len(nrow(top_genes))) {
    report <- c(report, sprintf("  %-15s : %d sites",
                                top_genes$gene_id[i], top_genes$N[i]))
  }

  report <- c(report, "")
  writeLines(report, con = output_file)
  message("Wrote report: ", output_file)

  invisible(report)
}
