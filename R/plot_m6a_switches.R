#' Plot All m6Aswitch Results — 5 Key Publication Plots
#'
#' Generates and saves the five core publication-quality figures from an
#' m6Aswitch analysis in a single call:
#'
#' \enumerate{
#'   \item \strong{m6A Fate Distribution} — bar chart of LOST / GAINED / RETAINED counts.
#'   \item \strong{LOST Mechanism Breakdown} — why m6A is lost: exon skipped vs unmethylated.
#'   \item \strong{Top Genes Heatmap} — top 30 genes by total m6A switching events.
#'   \item \strong{Volcano} — isoform switch magnitude (dIF) vs Δm6A probability.
#'   \item \strong{Individual Gene Track Plots} — per-gene m6A site tracks for the top 5 genes.
#' }
#'
#' @param classified data.table. Output of \code{classify_lost_mechanism()}.
#'   Must contain columns: \code{gene_id}, \code{isoform_a}, \code{isoform_b},
#'   \code{m6a_fate}, \code{lost_mechanism}, \code{probability_a},
#'   \code{probability_b}, \code{m6a_in_isoform_a}, \code{m6a_in_isoform_b},
#'   and \code{transcript_position}.
#' @param switch_pairs data.table. The isoform switch table used as input to
#'   the analysis (e.g. from \code{parse_isoform_switch()}). Must contain
#'   columns \code{gene_id}, \code{dif} (direction of isoform fraction change),
#'   and \code{fdr}. Used for the volcano plot (Plot 4) and gene-level subtitles
#'   (Plot 5). Pass \code{NULL} to skip those elements.
#' @param output_dir Character. Directory in which to write PDF files. Created
#'   recursively if it does not exist.
#' @param analysis_name Character. Label prepended to each plot title, e.g.
#'   \code{"Analysis 1: Astrocyte vs EGFRv3+IDHwt"} (default: \code{"Analysis"}).
#' @param top_n_genes Integer. Number of genes to show in the heatmap (default: 30).
#' @param top_n_gene_plots Integer. Number of top genes for which individual
#'   track plots (Plot 5) are generated (default: 5).
#'
#' @return Invisibly returns a named list of the ggplot objects that were
#'   successfully created:
#'   \code{p1_fate}, \code{p2_mechanism}, \code{p3_heatmap}, \code{p4_volcano},
#'   and \code{p5_genes} (a named list of per-gene plots).
#'
#' @examples
#' \dontrun{
#' # After running the full pipeline:
#' # Analysis 1
#' plot_m6aswitch_results(
#'   classified    = classified_a1,
#'   switch_pairs  = iso_switches_a1,
#'   output_dir    = "results/Analysis1/plots",
#'   analysis_name = "Analysis 1: Astrocyte vs EGFRv3+IDHwt"
#' )
#'
#' # Analysis 2
#' plot_m6aswitch_results(
#'   classified    = classified_a2,
#'   switch_pairs  = iso_switches_a2,
#'   output_dir    = "results/Analysis2/plots",
#'   analysis_name = "Analysis 2: Astrocyte vs All Tumors"
#' )
#' }
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_m6aswitch_results <- function(classified,
                                   switch_pairs  = NULL,
                                   output_dir,
                                   analysis_name   = "Analysis",
                                   top_n_genes     = 30L,
                                   top_n_gene_plots = 5L) {

  if (!data.table::is.data.table(classified) || nrow(classified) == 0) {
    stop("classified must be a non-empty data.table from classify_lost_mechanism()")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  dt <- data.table::copy(classified)

  fate_levels <- c("LOST", "GAINED", "RETAINED")
  fate_colors <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")

  dt[, m6a_fate := factor(m6a_fate, levels = fate_levels)]

  base_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40")
    )

  plots <- list()

  # ── Plot 1: m6A Fate Distribution ─────────────────────────────────────────
  message("Plot 1: m6A fate distribution...")

  fate_counts <- dt[, .N, by = m6a_fate][order(m6a_fate)]
  fate_counts[, pct := 100 * N / sum(N)]

  p1 <- ggplot2::ggplot(
    fate_counts,
    ggplot2::aes(x = m6a_fate, y = N, fill = m6a_fate)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
      vjust = -0.4, size = 4
    ) +
    ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
    ggplot2::labs(
      title    = paste(analysis_name, "\u2014 m6A Site Fate Across Isoform Switches"),
      subtitle = "What happens to m6A sites when genes switch isoforms?",
      x        = "m6A Fate",
      y        = "Number of m6A sites",
      fill     = "Fate"
    ) +
    ggplot2::expand_limits(y = max(fate_counts$N) * 1.15) +
    base_theme

  ggplot2::ggsave(
    file.path(output_dir, "plot1_m6a_fate_distribution.pdf"),
    p1, width = 8, height = 6
  )
  plots$p1_fate <- p1
  message("  Saved: plot1_m6a_fate_distribution.pdf")

  # ── Plot 2: LOST Mechanism Breakdown ──────────────────────────────────────
  message("Plot 2: LOST mechanism breakdown...")

  lost_dt <- dt[m6a_fate == "LOST" & !is.na(lost_mechanism)]

  if (nrow(lost_dt) == 0) {
    message("  WARNING: No LOST sites with mechanism assigned — skipping Plot 2.")
    message("  Run classify_lost_mechanism() first.")
  } else {
    mech_counts <- lost_dt[, .N, by = lost_mechanism][order(-N)]
    mech_counts[, pct         := 100 * N / sum(N)]
    mech_counts[, lost_mechanism := factor(lost_mechanism, levels = lost_mechanism)]

    mech_colors <- c(
      LOST_EXON_SKIPPED = "#E69F00",
      LOST_UNMETHYLATED = "#CC79A7"
    )

    p2 <- ggplot2::ggplot(
      mech_counts,
      ggplot2::aes(x = lost_mechanism, y = N, fill = lost_mechanism)
    ) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
        vjust = -0.4, size = 4
      ) +
      ggplot2::scale_fill_manual(values = mech_colors, drop = FALSE) +
      ggplot2::labs(
        title    = paste(analysis_name, "\u2014 Why is m6A LOST During Isoform Switching?"),
        subtitle = paste0(
          "EXON_SKIPPED: exon carrying m6A is absent from tumor isoform\n",
          "UNMETHYLATED: exon present in tumor isoform but no longer methylated"
        ),
        x        = "Mechanism",
        y        = "Number of LOST m6A sites",
        fill     = "Mechanism"
      ) +
      ggplot2::expand_limits(y = max(mech_counts$N) * 1.15) +
      base_theme +
      ggplot2::theme(legend.position = "none")

    ggplot2::ggsave(
      file.path(output_dir, "plot2_lost_mechanism.pdf"),
      p2, width = 8, height = 6
    )
    plots$p2_mechanism <- p2
    message("  Saved: plot2_lost_mechanism.pdf")
  }

  # ── Plot 3: Top Genes Heatmap ──────────────────────────────────────────────
  message("Plot 3: Top genes heatmap...")

  gene_totals <- dt[, .N, by = gene_id][order(-N)]
  keep_genes  <- head(gene_totals$gene_id, top_n_genes)

  gene_fate <- dt[gene_id %in% keep_genes, .N, by = .(gene_id, m6a_fate)]
  gene_fate[, gene_id := factor(gene_id, levels = rev(keep_genes))]

  p3 <- ggplot2::ggplot(
    gene_fate,
    ggplot2::aes(x = m6a_fate, y = gene_id, fill = N)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = N), size = 3, color = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_viridis_c(option = "C", name = "Count") +
    ggplot2::labs(
      title    = paste(analysis_name, sprintf(
        "\u2014 Top %d Genes by m6A Switching Events", length(keep_genes)
      )),
      subtitle = "Number of m6A sites per fate category per gene",
      x        = "m6A Fate",
      y        = "Gene"
    ) +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

  ggplot2::ggsave(
    file.path(output_dir, "plot3_top_genes_heatmap.pdf"),
    p3, width = 9, height = max(7, 0.32 * length(keep_genes) + 2)
  )
  plots$p3_heatmap <- p3
  message("  Saved: plot3_top_genes_heatmap.pdf")

  # ── Plot 4: Volcano — dIF vs Δm6A Probability ─────────────────────────────
  message("Plot 4: Volcano dIF vs m6A change...")

  if (is.null(switch_pairs)) {
    message("  WARNING: switch_pairs is NULL — skipping Plot 4.")
  } else {
    if (!data.table::is.data.table(switch_pairs)) {
      switch_pairs <- data.table::as.data.table(switch_pairs)
    }

    # Flexible column lookup: accept gene_id or geneID
    sp <- data.table::copy(switch_pairs)
    if ("geneID" %in% names(sp) && !"gene_id" %in% names(sp)) {
      data.table::setnames(sp, "geneID", "gene_id")
    }
    if ("dIF" %in% names(sp) && !"dif" %in% names(sp)) {
      data.table::setnames(sp, "dIF", "dif")
    }
    if ("isoform_switch_q_value" %in% names(sp) && !"fdr" %in% names(sp)) {
      data.table::setnames(sp, "isoform_switch_q_value", "fdr")
    }

    volcano_dt <- data.table::copy(dt)
    if (!"dif" %in% names(volcano_dt)) {
      merge_cols <- intersect(c("gene_id", "dif"), names(sp))
      if (all(c("gene_id", "dif") %in% merge_cols)) {
        dif_lookup <- unique(sp[, ..merge_cols], by = "gene_id")
        volcano_dt <- merge(
          volcano_dt,
          dif_lookup,
          by = "gene_id",
          all.x = TRUE,
          sort = FALSE
        )
      }
    }

    volcano_dt[, delta_prob := data.table::fcoalesce(probability_b, 0) -
                             data.table::fcoalesce(probability_a, 0)]

    plot_data <- volcano_dt[!is.na(delta_prob) & !is.na(dif)]

    if (nrow(plot_data) == 0) {
      message("  WARNING: No rows with both delta_prob and dif — skipping Plot 4.")
    } else {
      p4 <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = dif, y = delta_prob, color = m6a_fate)
      ) +
        ggplot2::geom_point(alpha = 0.45, size = 1.5) +
        ggplot2::geom_hline(yintercept = 0,          linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = 0,          linetype = "dashed",
                            linewidth = 0.4) +
        ggplot2::geom_vline(xintercept = c(-0.1, 0.1),
                            linetype = "dotted", linewidth = 0.3, color = "grey60") +
        ggplot2::scale_color_manual(values = fate_colors, drop = FALSE) +
        ggplot2::labs(
          title    = paste(analysis_name,
                           "\u2014 Isoform Switch Magnitude vs m6A Change"),
          subtitle = paste0(
            "X: how much isoform usage changes between conditions | ",
            "Y: how much m6A probability changes (not detected treated as 0)"
          ),
          x        = "dIF (isoform fraction difference, tumor \u2212 normal)",
          y        = expression(Delta * " m6A probability (isoform B \u2212 isoform A)"),
          color    = "m6A Fate"
        ) +
        base_theme

      ggplot2::ggsave(
        file.path(output_dir, "plot4_volcano_dIF_vs_m6A.pdf"),
        p4, width = 9, height = 7
      )
      plots$p4_volcano <- p4
      message("  Saved: plot4_volcano_dIF_vs_m6A.pdf")
    }
  }

  # ── Plot 5: Individual Gene Track Plots ───────────────────────────────────
  message("Plot 5: Individual gene track plots (top ", top_n_gene_plots, " genes)...")

  top_genes <- head(gene_totals$gene_id, top_n_gene_plots)
  plots$p5_genes <- list()

  for (gene in top_genes) {

    gene_dt <- dt[gene_id == gene]
    if (nrow(gene_dt) == 0) next

    # Determine x variable: prefer genomic coordinate for comparable isoform tracks
    x_col <- if ("start" %in% names(gene_dt)) {
      "start"
    } else if ("transcript_position" %in% names(gene_dt)) {
      "transcript_position"
    } else if ("position" %in% names(gene_dt)) {
      "position"
    } else {
      message("  Skipping ", gene, ": no position column found.")
      next
    }

    dt_a <- gene_dt[, .(
      isoform  = isoform_a,
      position = get(x_col),
      m6a_fate = m6a_fate,
      present  = m6a_in_isoform_a,
      side     = "Isoform A (Normal / Astrocyte)"
    )]
    dt_b <- gene_dt[, .(
      isoform  = isoform_b,
      position = get(x_col),
      m6a_fate = m6a_fate,
      present  = m6a_in_isoform_b,
      side     = "Isoform B (Tumor)"
    )]

    long_dt <- data.table::rbindlist(list(dt_a, dt_b), use.names = TRUE)
    long_dt  <- long_dt[present == TRUE]

    if (nrow(long_dt) == 0) {
      message("  Skipping ", gene, ": no m6A-positive positions to plot.")
      next
    }

    # Build subtitle — pull dIF and q-value from switch_pairs if available
    subtitle_text <- sprintf("Switch: Astrocyte \u2192 Tumor | Top gene by m6A switching events")
    if (!is.null(switch_pairs)) {
      sp2 <- data.table::copy(switch_pairs)
      if ("geneID"                %in% names(sp2)) data.table::setnames(sp2, "geneID",                "gene_id", skip_absent = TRUE)
      if ("dIF"                   %in% names(sp2)) data.table::setnames(sp2, "dIF",                   "dif",     skip_absent = TRUE)
      if ("isoform_switch_q_value"%in% names(sp2)) data.table::setnames(sp2, "isoform_switch_q_value","fdr",     skip_absent = TRUE)

      gene_row <- sp2[gene_id == gene]
      if (nrow(gene_row) > 0) {
        dif_val  <- gene_row$dif[1]
        fdr_val  <- gene_row$fdr[1]
        subtitle_text <- sprintf(
          "dIF = %.2f | FDR = %.2e | Switch: Astrocyte \u2192 Tumor",
          dif_val, fdr_val
        )
      }
    }

    x_range_pad <- diff(range(long_dt$position, na.rm = TRUE)) * 0.05
    x_min <- min(long_dt$position, na.rm = TRUE) - x_range_pad
    x_max <- max(long_dt$position, na.rm = TRUE) + x_range_pad

    p5 <- ggplot2::ggplot(
      long_dt,
      ggplot2::aes(x = position, y = side, color = m6a_fate)
    ) +
      # Backbone representing the isoform (grey bar)
      ggplot2::geom_segment(
        ggplot2::aes(x = x_min, xend = x_max, y = side, yend = side),
        color = "grey80", linewidth = 3, inherit.aes = FALSE,
        data = unique(long_dt[, .(side)])
      ) +
      # m6A site points
      ggplot2::geom_point(size = 4, alpha = 0.95) +
      # Fate label above each point
      ggplot2::geom_text(
        ggplot2::aes(label = m6a_fate),
        vjust = -1.2, size = 2.8, fontface = "bold"
      ) +
      ggplot2::scale_color_manual(values = fate_colors, drop = FALSE) +
      ggplot2::labs(
        title    = paste0(gene, " \u2014 m6A Sites on Switching Isoforms"),
        subtitle = subtitle_text,
        x        = if (x_col == "start") "Genomic coordinate" else if (x_col == "transcript_position") "Transcript position (nt)" else "Position",
        y        = "",
        color    = "m6A Fate"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        legend.position  = "right",
        axis.text.y      = ggplot2::element_text(size = 11, face = "bold"),
        plot.title       = ggplot2::element_text(face = "bold", size = 13),
        plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40"),
        panel.border     = ggplot2::element_rect(fill = NA, color = "grey85")
      )

    fname <- file.path(
      output_dir,
      paste0("plot5_gene_", gsub("[^A-Za-z0-9_.-]", "_", gene), ".pdf")
    )
    ggplot2::ggsave(fname, p5, width = 10, height = 5)
    plots$p5_genes[[gene]] <- p5
    message("  Saved: plot5_gene_", gene, ".pdf")
  }

  # ── Summary ────────────────────────────────────────────────────────────────
  message("\n=== All plots saved to: ", output_dir, " ===")
  message("Files:")
  for (f in list.files(output_dir, pattern = "\\.pdf$", full.names = FALSE)) {
    message("  ", f)
  }

  invisible(plots)
}


# ── Legacy dispatcher (kept for backwards compatibility) ────────────────────

#' Plot m6A Switching Events (legacy dispatcher)
#'
#' @description
#' Superseded.
#'
#' This function is superseded by \code{\link{plot_m6aswitch_results}()}, which
#' generates all five publication plots in one call and saves them to disk.
#'
#' \code{plot_m6a_switches()} is kept for backwards compatibility and exploratory
#' use; it returns a single ggplot object for the specified \code{plot_type}.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#' @param plot_type One of \code{"summary"}, \code{"heatmap"}, \code{"delta_prob"}.
#' @param top_n_genes Integer. Number of genes to display in heatmap (default: 30).
#'
#' @return ggplot object.
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_m6a_switches <- function(m6a_switches,
                              plot_type   = "summary",
                              top_n_genes = 30) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  plot_type <- match.arg(plot_type, c("summary", "heatmap", "delta_prob"))

  fate_levels <- c("LOST", "GAINED", "RETAINED")
  fate_colors <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")

  dt <- data.table::copy(m6a_switches)
  dt[, m6a_fate := factor(m6a_fate, levels = fate_levels)]

  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right"
    )

  if (plot_type == "summary") {
    fate_counts <- dt[, .N, by = m6a_fate][order(m6a_fate)]
    fate_counts[, pct := 100 * N / sum(N)]

    p <- ggplot2::ggplot(
      fate_counts,
      ggplot2::aes(x = m6a_fate, y = N, fill = m6a_fate)
    ) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
        vjust = -0.4, size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
      ggplot2::labs(
        title = "m6A Site Fate Across Isoform Switches",
        x = "m6A fate", y = "Number of sites", fill = "Fate"
      ) +
      ggplot2::expand_limits(y = max(fate_counts$N) * 1.12) +
      base_theme

    return(p)
  }

  if (plot_type == "heatmap") {
    gene_totals <- dt[, .N, by = gene_id][order(-N)]
    keep_genes  <- head(gene_totals$gene_id, top_n_genes)

    gene_fate <- dt[gene_id %in% keep_genes, .N, by = .(gene_id, m6a_fate)]
    gene_fate[, gene_id := factor(gene_id, levels = rev(keep_genes))]

    p <- ggplot2::ggplot(
      gene_fate,
      ggplot2::aes(x = m6a_fate, y = gene_id, fill = N)
    ) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Count") +
      ggplot2::labs(
        title = sprintf("m6A Fate by Gene (Top %d genes)", length(keep_genes)),
        x = "m6A fate", y = "Gene"
      ) +
      base_theme

    return(p)
  }

  if (plot_type == "delta_prob") {
    dt[, delta_prob := data.table::fcoalesce(probability_b, 0) -
                        data.table::fcoalesce(probability_a, 0)]

    p <- ggplot2::ggplot(
      dt[!is.na(delta_prob)],
      ggplot2::aes(x = m6a_fate, y = delta_prob, fill = m6a_fate)
    ) +
      ggplot2::geom_violin(alpha = 0.6, trim = TRUE, color = NA) +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.9) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
      ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
      ggplot2::labs(
        title = "Change in m6A Probability Across Switches",
        subtitle = "Not-detected sites are treated as probability 0 for delta calculations",
        x = "m6A fate",
        y = expression(Delta * " probability (isoform B - isoform A)"),
        fill = "Fate"
      ) +
      base_theme

    return(p)
  }
}


#' Plot Isoform Switch Details for a Single Gene
#'
#' Creates a lollipop/track-style plot showing m6A sites on both switching
#' isoforms for one gene. Points are coloured by m6A fate (LOST / GAINED /
#' RETAINED) and positioned by genomic coordinate when available.
#'
#' @param gene_id Character. Gene name or ID to plot.
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#' @param iso_lengths Optional named numeric vector: \code{isoform_id -> length}.
#'   When supplied, positions are expressed as a fraction [0, 1] of total
#'   isoform length.
#' @param switch_pairs Optional data.table from \code{parse_isoform_switch()}.
#'   When supplied, dIF and FDR are added to the plot subtitle.
#'
#' @return ggplot object.
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_isoform_details <- function(gene_id,
                                 m6a_switches,
                                 iso_lengths  = NULL,
                                 switch_pairs = NULL) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  # Avoid data.table NSE conflict with function argument name
  .gene <- gene_id
  gene_data <- m6a_switches[gene_id == .gene]

  if (nrow(gene_data) == 0) {
    stop(sprintf("No data for gene: %s", gene_id))
  }

  fate_colors <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")

  # Determine position column: prefer genomic coordinate for cross-isoform comparison
  x_col <- if ("start" %in% names(gene_data)) {
    "start"
  } else if ("transcript_position" %in% names(gene_data)) {
    "transcript_position"
  } else {
    "position"
  }

  dt_a <- gene_data[, .(
    isoform  = isoform_a,
    position = get(x_col),
    m6a_fate = m6a_fate,
    present  = m6a_in_isoform_a,
    side     = "Isoform A (Normal / Astrocyte)"
  )]
  dt_b <- gene_data[, .(
    isoform  = isoform_b,
    position = get(x_col),
    m6a_fate = m6a_fate,
    present  = m6a_in_isoform_b,
    side     = "Isoform B (Tumor)"
  )]

  long_dt <- data.table::rbindlist(list(dt_a, dt_b), use.names = TRUE)
  long_dt  <- long_dt[present == TRUE]

  if (nrow(long_dt) == 0) {
    stop(sprintf("No m6A-present positions to plot for gene: %s", gene_id))
  }

  # Optional scaling by isoform length
  if (!is.null(iso_lengths)) {
    if (is.null(names(iso_lengths))) {
      stop("iso_lengths must be a named vector: names are isoform IDs")
    }
    long_dt[, iso_length := as.numeric(iso_lengths[isoform])]
    long_dt[, position   := ifelse(
      !is.na(iso_length) & iso_length > 0,
      position / iso_length,
      NA_real_
    )]
    x_lab <- "Relative position within isoform (0 = 5', 1 = 3')"
  } else {
    x_lab <- if (x_col == "start") "Genomic coordinate" else if (x_col == "transcript_position") "Transcript position (nt)" else "Position"
  }

  subtitle_text <- sprintf("Switch: Astrocyte \u2192 Tumor")
  if (!is.null(switch_pairs)) {
    sp <- data.table::copy(switch_pairs)
    if ("geneID"                %in% names(sp)) data.table::setnames(sp, "geneID",                "gene_id", skip_absent = TRUE)
    if ("dIF"                   %in% names(sp)) data.table::setnames(sp, "dIF",                   "dif",     skip_absent = TRUE)
    if ("isoform_switch_q_value"%in% names(sp)) data.table::setnames(sp, "isoform_switch_q_value","fdr",     skip_absent = TRUE)
    row <- sp[gene_id == .gene]
    if (nrow(row) > 0) {
      subtitle_text <- sprintf(
        "dIF = %.2f | FDR = %.2e | Switch: Astrocyte \u2192 Tumor",
        row$dif[1], row$fdr[1]
      )
    }
  }

  p <- ggplot2::ggplot(
    long_dt,
    ggplot2::aes(x = position, y = side, color = m6a_fate)
  ) +
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    ggplot2::scale_color_manual(values = fate_colors, drop = FALSE) +
    ggplot2::labs(
      title    = sprintf("m6A Sites in Isoform Switches: %s", gene_id),
      subtitle = subtitle_text,
      x        = x_lab,
      y        = "",
      color    = "Fate"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      axis.text.y      = ggplot2::element_text(size = 11, face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold")
    )

  return(p)
}
