#' Plot All m6Aswitch Results - 5 Key Publication Plots
#'
#' Generates and saves the five core publication figures from an m6Aswitch
#' analysis in a single call:
#'
#' \enumerate{
#'   \item \strong{Isoform Status Distribution} - bar chart of
#'     ISOFORM_A_ONLY / ISOFORM_B_ONLY / IN_BOTH_ISOFORMS counts.
#'   \item \strong{Mechanism Breakdown} - structural (exon skipped) vs
#'     regulatory (present but unmethylated).
#'   \item \strong{Top Genes Heatmap} - top N genes by total m6A sites.
#'   \item \strong{Volcano} - isoform switch magnitude (dIF) vs delta m6A
#'     probability.
#'   \item \strong{Individual Gene Track Plots} - per-gene m6A site tracks.
#' }
#'
#' @param classified data.table from \code{classify_lost_mechanism()}.
#'   Must contain: \code{gene_id}, \code{isoform_a}, \code{isoform_b},
#'   \code{isoform_status}, \code{isoform_mechanism}, \code{probability_a},
#'   \code{probability_b}, \code{m6a_in_isoform_a}, \code{m6a_in_isoform_b},
#'   and a position column.
#' @param switch_pairs data.table from \code{parse_isoform_switch()}, used for
#'   the volcano plot and gene subtitles. Pass \code{NULL} to skip those.
#' @param output_dir Directory for PDF output. Created if absent.
#' @param analysis_name Label prepended to each plot title.
#' @param top_n_genes Number of genes in the heatmap (default 30).
#' @param top_n_gene_plots Number of per-gene track plots (default 5).
#' @param color_by Which classification to colour by: \code{"isoform_status"}
#'   (structural, the default) or \code{"m6a_fate"} (regulatory - requires that
#'   condition information was supplied upstream).
#'
#' @return Invisibly, a named list of the ggplot objects created.
#'
#' @section Structural vs regulatory:
#' \code{isoform_status} describes which transcript of the switch pair carries
#' each site. \code{m6a_fate} describes which condition detected it, and is NA
#' unless a \code{condition} column was supplied to
#' \code{lift_m6a_to_genomic()}. These are different questions - see
#' \code{\link{annotate_m6a_switches_genomic}}.
#'
#' @examples
#' \dontrun{
#' plot_m6aswitch_results(
#'   classified    = classified,
#'   switch_pairs  = iso_switches,
#'   output_dir    = "results/plots",
#'   analysis_name = "Astrocyte vs Tumor"
#' )
#'
#' # colour by condition instead (requires condition column upstream)
#' plot_m6aswitch_results(classified, iso_switches, "results/plots_cond",
#'                        color_by = "m6a_fate")
#' }
#'
#' @import ggplot2
#' @import data.table
#' @importFrom utils head
#'
#' @export
plot_m6aswitch_results <- function(classified,
                                   switch_pairs     = NULL,
                                   output_dir,
                                   analysis_name    = "Analysis",
                                   top_n_genes      = 30L,
                                   top_n_gene_plots = 5L,
                                   color_by         = c("isoform_status", "m6a_fate")) {

  if (!data.table::is.data.table(classified) || nrow(classified) == 0) {
    stop("classified must be a non-empty data.table from classify_lost_mechanism()")
  }

  color_by <- match.arg(color_by)

  if (!color_by %in% names(classified)) {
    stop("Column '", color_by, "' not found in classified. Available: ",
         paste(names(classified), collapse = ", "))
  }

  if (color_by == "m6a_fate" && all(is.na(classified$m6a_fate))) {
    stop("m6a_fate is entirely NA. Condition information was not supplied to ",
         "lift_m6a_to_genomic(). Use color_by = 'isoform_status' instead.")
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  dt <- data.table::copy(classified)

  # ── Level sets and colours for whichever classification is being plotted ────
  if (color_by == "isoform_status") {
    lvls <- c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY", "IN_BOTH_ISOFORMS")
    cols <- c(ISOFORM_A_ONLY   = "#D55E00",
              ISOFORM_B_ONLY   = "#009E73",
              IN_BOTH_ISOFORMS = "#0072B2")
    legend_lab   <- "Isoform status"
    axis_lab     <- "Isoform status"
    plot1_sub    <- "Which isoform of the switch pair carries each m6A site?"
  } else {
    lvls <- c("LOST", "GAINED", "RETAINED")
    cols <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")
    legend_lab   <- "m6A fate"
    axis_lab     <- "m6A fate"
    plot1_sub    <- "Which condition detected each m6A site?"
  }

  dt[, .plot_class := factor(get(color_by), levels = lvls)]
  dt <- dt[!is.na(.plot_class)]

  if (nrow(dt) == 0) {
    stop("No rows with a non-NA value in '", color_by, "'.")
  }

  base_theme <- ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      plot.title       = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle    = ggplot2::element_text(size = 10, color = "grey40")
    )

  plots <- list()

  # ── Plot 1: distribution ────────────────────────────────────────────────────
  message("Plot 1: ", color_by, " distribution...")

  fate_counts <- dt[, .N, by = .plot_class][order(.plot_class)]
  fate_counts[, pct := 100 * N / sum(N)]

  p1 <- ggplot2::ggplot(fate_counts,
                        ggplot2::aes(x = .plot_class, y = N, fill = .plot_class)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
                       vjust = -0.4, size = 4) +
    ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
    ggplot2::labs(
      title    = paste(analysis_name, "\u2014 m6A Sites Across Isoform Switches"),
      subtitle = plot1_sub,
      x = axis_lab, y = "Number of m6A sites", fill = legend_lab
    ) +
    ggplot2::expand_limits(y = max(fate_counts$N) * 1.15) +
    base_theme

  ggplot2::ggsave(file.path(output_dir, "plot1_distribution.pdf"),
                  p1, width = 8, height = 6)
  plots$p1_fate <- p1
  message("  Saved: plot1_distribution.pdf")

  # ── Plot 2: mechanism breakdown ─────────────────────────────────────────────
  message("Plot 2: mechanism breakdown...")

  if (!"isoform_mechanism" %in% names(dt)) {
    message("  No isoform_mechanism column \u2014 skipping Plot 2. ",
            "Run classify_lost_mechanism() first.")
  } else {
    mech_dt <- dt[!is.na(isoform_mechanism)]

    if (nrow(mech_dt) == 0) {
      message("  No sites with a mechanism assigned \u2014 skipping Plot 2.")
    } else {
      mech_counts <- mech_dt[, .N, by = isoform_mechanism][order(-N)]
      mech_counts[, pct := 100 * N / sum(N)]
      mech_counts[, isoform_mechanism := factor(isoform_mechanism,
                                                levels = isoform_mechanism)]

      mech_colors <- c(
        A_ONLY_EXON_SKIPPED    = "#E69F00",
        A_ONLY_UNMETHYLATED    = "#CC79A7",
        B_ONLY_EXON_INCLUDED   = "#56B4E9",
        B_ONLY_NEW_METHYLATION = "#009E73"
      )

      p2 <- ggplot2::ggplot(mech_counts,
                            ggplot2::aes(x = isoform_mechanism, y = N,
                                         fill = isoform_mechanism)) +
        ggplot2::geom_col(width = 0.6) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
                           vjust = -0.4, size = 4) +
        ggplot2::scale_fill_manual(values = mech_colors, drop = FALSE) +
        ggplot2::labs(
          title    = paste(analysis_name,
                           "\u2014 Why Do m6A Sites Differ Between Isoforms?"),
          subtitle = paste0(
            "EXON_SKIPPED / EXON_INCLUDED: structural \u2014 the sequence is absent from one isoform\n",
            "UNMETHYLATED / NEW_METHYLATION: regulatory \u2014 sequence present in both, modification on one"
          ),
          x = "Mechanism", y = "Number of sites", fill = "Mechanism"
        ) +
        ggplot2::expand_limits(y = max(mech_counts$N) * 1.15) +
        base_theme +
        ggplot2::theme(legend.position = "none",
                       axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

      ggplot2::ggsave(file.path(output_dir, "plot2_mechanism.pdf"),
                      p2, width = 9, height = 6)
      plots$p2_mechanism <- p2
      message("  Saved: plot2_mechanism.pdf")
    }
  }

  # ── Plot 3: top genes heatmap ───────────────────────────────────────────────
  message("Plot 3: top genes heatmap...")

  gene_totals <- dt[, .N, by = gene_id][order(-N)]
  keep_genes  <- utils::head(gene_totals$gene_id, top_n_genes)

  gene_fate <- dt[gene_id %in% keep_genes, .N, by = .(gene_id, .plot_class)]
  gene_fate[, gene_id := factor(gene_id, levels = rev(keep_genes))]

  p3 <- ggplot2::ggplot(gene_fate,
                        ggplot2::aes(x = .plot_class, y = gene_id, fill = N)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = N), size = 3,
                       color = "white", fontface = "bold") +
    ggplot2::scale_fill_viridis_c(option = "C", name = "Count") +
    ggplot2::labs(
      title    = paste(analysis_name,
                       sprintf("\u2014 Top %d Genes by m6A Sites", length(keep_genes))),
      subtitle = paste("Number of m6A sites per", legend_lab, "per gene"),
      x = axis_lab, y = "Gene"
    ) +
    base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9),
                   axis.text.x = ggplot2::element_text(angle = 20, hjust = 1))

  ggplot2::ggsave(file.path(output_dir, "plot3_top_genes_heatmap.pdf"),
                  p3, width = 9, height = max(7, 0.32 * length(keep_genes) + 2))
  plots$p3_heatmap <- p3
  message("  Saved: plot3_top_genes_heatmap.pdf")

  # ── Plot 4: volcano ─────────────────────────────────────────────────────────
  message("Plot 4: volcano dIF vs m6A change...")

  if (is.null(switch_pairs)) {
    message("  switch_pairs is NULL \u2014 skipping Plot 4.")
  } else {
    if (!data.table::is.data.table(switch_pairs)) {
      switch_pairs <- data.table::as.data.table(switch_pairs)
    }

    sp <- data.table::copy(switch_pairs)
    if ("geneID" %in% names(sp) && !"gene_id" %in% names(sp)) {
      data.table::setnames(sp, "geneID", "gene_id")
    }
    if ("isoformID_A" %in% names(sp) && !"isoform_a" %in% names(sp)) {
      data.table::setnames(sp, "isoformID_A", "isoform_a")
    }
    if ("isoformID_B" %in% names(sp) && !"isoform_b" %in% names(sp)) {
      data.table::setnames(sp, "isoformID_B", "isoform_b")
    }
    if ("dIF" %in% names(sp) && !"dif" %in% names(sp)) {
      data.table::setnames(sp, "dIF", "dif")
    }
    if ("isoform_switch_q_value" %in% names(sp) && !"fdr" %in% names(sp)) {
      data.table::setnames(sp, "isoform_switch_q_value", "fdr")
    }

    volcano_dt <- data.table::copy(dt)

    # Merge on the full switch pair, not gene alone - a gene with several
    # switch pairs would otherwise get one pair's dIF applied to all its sites.
    if (!"dif" %in% names(volcano_dt)) {
      key_cols <- intersect(c("gene_id", "isoform_a", "isoform_b"), names(sp))
      if (length(key_cols) == 3 && "dif" %in% names(sp)) {
        volcano_dt <- merge(volcano_dt,
                            unique(sp[, c(key_cols, "dif"), with = FALSE]),
                            by = key_cols, all.x = TRUE, sort = FALSE)
      } else if ("dif" %in% names(sp)) {
        message("  Falling back to gene-level dIF merge \u2014 ",
                "isoform columns not found in switch_pairs.")
        volcano_dt <- merge(volcano_dt,
                            unique(sp[, .(gene_id, dif)], by = "gene_id"),
                            by = "gene_id", all.x = TRUE, sort = FALSE)
      }
    }

    if (!"dif" %in% names(volcano_dt)) {
      message("  No dIF column available \u2014 skipping Plot 4.")
    } else {
      volcano_dt[, delta_prob := data.table::fcoalesce(probability_b, 0) -
                                 data.table::fcoalesce(probability_a, 0)]
      plot_data <- volcano_dt[!is.na(delta_prob) & !is.na(dif)]

      if (nrow(plot_data) == 0) {
        message("  No rows with both delta_prob and dif \u2014 skipping Plot 4.")
      } else {
        p4 <- ggplot2::ggplot(plot_data,
                              ggplot2::aes(x = dif, y = delta_prob,
                                           color = .plot_class)) +
          ggplot2::geom_point(alpha = 0.45, size = 1.5) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
          ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
          ggplot2::geom_vline(xintercept = c(-0.1, 0.1), linetype = "dotted",
                              linewidth = 0.3, color = "grey60") +
          ggplot2::scale_color_manual(values = cols, drop = FALSE) +
          ggplot2::labs(
            title    = paste(analysis_name,
                             "\u2014 Isoform Switch Magnitude vs m6A Difference"),
            subtitle = paste0("X: change in isoform usage between conditions | ",
                              "Y: probability difference between isoforms ",
                              "(undetected treated as 0)"),
            x     = "dIF (isoform fraction difference)",
            y     = expression(Delta * " m6A probability (isoform B \u2212 isoform A)"),
            color = legend_lab
          ) +
          base_theme

        ggplot2::ggsave(file.path(output_dir, "plot4_volcano_dIF_vs_m6A.pdf"),
                        p4, width = 9, height = 7)
        plots$p4_volcano <- p4
        message("  Saved: plot4_volcano_dIF_vs_m6A.pdf")
      }
    }
  }

  # ── Plot 5: per-gene tracks ─────────────────────────────────────────────────
  message("Plot 5: per-gene track plots (top ", top_n_gene_plots, ")...")

  top_genes      <- utils::head(gene_totals$gene_id, top_n_gene_plots)
  plots$p5_genes <- list()

  for (gene in top_genes) {

    gene_dt <- dt[gene_id == gene]
    if (nrow(gene_dt) == 0) next

    x_col <- if ("start" %in% names(gene_dt)) {
      "start"
    } else if ("transcript_position" %in% names(gene_dt)) {
      "transcript_position"
    } else if ("position" %in% names(gene_dt)) {
      "position"
    } else {
      message("  Skipping ", gene, ": no position column.")
      next
    }

    cond_a_lab <- if (!is.na(gene_dt$condition_1[1])) gene_dt$condition_1[1] else "condition 1"
    cond_b_lab <- if (!is.na(gene_dt$condition_2[1])) gene_dt$condition_2[1] else "condition 2"

    dt_a <- gene_dt[, .(
      isoform      = isoform_a,
      position     = get(x_col),
      .plot_class  = .plot_class,
      present      = m6a_in_isoform_a,
      side         = paste0("Isoform A (", cond_a_lab, ")")
    )]
    dt_b <- gene_dt[, .(
      isoform      = isoform_b,
      position     = get(x_col),
      .plot_class  = .plot_class,
      present      = m6a_in_isoform_b,
      side         = paste0("Isoform B (", cond_b_lab, ")")
    )]

    long_dt <- data.table::rbindlist(list(dt_a, dt_b), use.names = TRUE)
    long_dt <- long_dt[present == TRUE]

    if (nrow(long_dt) == 0) {
      message("  Skipping ", gene, ": no m6A-positive positions.")
      next
    }

    subtitle_text <- sprintf("%s \u2192 %s", cond_a_lab, cond_b_lab)
    if (!is.na(gene_dt$dif[1]) && !is.na(gene_dt$fdr[1])) {
      subtitle_text <- sprintf("dIF = %.2f | FDR = %.2e | %s \u2192 %s",
                               gene_dt$dif[1], gene_dt$fdr[1],
                               cond_a_lab, cond_b_lab)
    }

    rng   <- range(long_dt$position, na.rm = TRUE)
    pad   <- max(diff(rng) * 0.05, 1)
    x_min <- rng[1] - pad
    x_max <- rng[2] + pad

    p5 <- ggplot2::ggplot(long_dt,
                          ggplot2::aes(x = position, y = side, color = .plot_class)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = x_min, xend = x_max, y = side, yend = side),
        color = "grey80", linewidth = 3, inherit.aes = FALSE,
        data = unique(long_dt[, .(side)])
      ) +
      ggplot2::geom_point(size = 4, alpha = 0.95) +
      ggplot2::scale_color_manual(values = cols, drop = FALSE) +
      ggplot2::labs(
        title    = paste0(gene, " \u2014 m6A Sites on Switching Isoforms"),
        subtitle = subtitle_text,
        x        = if (x_col == "start") "Genomic coordinate"
                   else if (x_col == "transcript_position") "Transcript position (nt)"
                   else "Position",
        y        = "",
        color    = legend_lab
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

    fname <- file.path(output_dir,
                       paste0("plot5_gene_",
                              gsub("[^A-Za-z0-9_.-]", "_", gene), ".pdf"))
    ggplot2::ggsave(fname, p5, width = 10, height = 5)
    plots$p5_genes[[gene]] <- p5
    message("  Saved: plot5_gene_", gene, ".pdf")
  }

  message("\n=== All plots saved to: ", output_dir, " ===")
  for (f in list.files(output_dir, pattern = "\\.pdf$")) message("  ", f)

  invisible(plots)
}


# ── Legacy dispatcher ────────────────────────────────────────────────────────

#' Plot m6A Switching Events (legacy dispatcher)
#'
#' @description
#' Superseded by \code{\link{plot_m6aswitch_results}()}, which generates all
#' five publication plots in one call. Kept for exploratory use; returns a
#' single ggplot object.
#'
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#' @param plot_type One of \code{"summary"}, \code{"heatmap"}, \code{"delta_prob"}.
#' @param top_n_genes Number of genes in the heatmap (default 30).
#' @param color_by \code{"isoform_status"} (default) or \code{"m6a_fate"}.
#'
#' @return ggplot object.
#'
#' @import ggplot2
#' @import data.table
#' @importFrom utils head
#'
#' @export
plot_m6a_switches <- function(m6a_switches,
                              plot_type   = "summary",
                              top_n_genes = 30,
                              color_by    = c("isoform_status", "m6a_fate")) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  plot_type <- match.arg(plot_type, c("summary", "heatmap", "delta_prob"))
  color_by  <- match.arg(color_by)

  if (!color_by %in% names(m6a_switches)) {
    stop("Column '", color_by, "' not found in m6a_switches.")
  }

  if (color_by == "isoform_status") {
    lvls <- c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY", "IN_BOTH_ISOFORMS")
    cols <- c(ISOFORM_A_ONLY   = "#D55E00",
              ISOFORM_B_ONLY   = "#009E73",
              IN_BOTH_ISOFORMS = "#0072B2")
    lab  <- "Isoform status"
  } else {
    lvls <- c("LOST", "GAINED", "RETAINED")
    cols <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")
    lab  <- "m6A fate"
  }

  dt <- data.table::copy(m6a_switches)
  dt[, .plot_class := factor(get(color_by), levels = lvls)]
  dt <- dt[!is.na(.plot_class)]

  if (nrow(dt) == 0) {
    stop("No rows with a non-NA value in '", color_by, "'.")
  }

  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.position  = "right")

  if (plot_type == "summary") {
    fate_counts <- dt[, .N, by = .plot_class][order(.plot_class)]
    fate_counts[, pct := 100 * N / sum(N)]

    return(
      ggplot2::ggplot(fate_counts,
                      ggplot2::aes(x = .plot_class, y = N, fill = .plot_class)) +
        ggplot2::geom_col(width = 0.75) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
                           vjust = -0.4, size = 3.5) +
        ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
        ggplot2::labs(title = "m6A Sites Across Isoform Switches",
                      x = lab, y = "Number of sites", fill = lab) +
        ggplot2::expand_limits(y = max(fate_counts$N) * 1.12) +
        base_theme
    )
  }

  if (plot_type == "heatmap") {
    gene_totals <- dt[, .N, by = gene_id][order(-N)]
    keep_genes  <- utils::head(gene_totals$gene_id, top_n_genes)

    gene_fate <- dt[gene_id %in% keep_genes, .N, by = .(gene_id, .plot_class)]
    gene_fate[, gene_id := factor(gene_id, levels = rev(keep_genes))]

    return(
      ggplot2::ggplot(gene_fate,
                      ggplot2::aes(x = .plot_class, y = gene_id, fill = N)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.2) +
        ggplot2::scale_fill_viridis_c(option = "C", name = "Count") +
        ggplot2::labs(
          title = sprintf("m6A Sites by Gene (Top %d)", length(keep_genes)),
          x = lab, y = "Gene") +
        base_theme
    )
  }

  if (plot_type == "delta_prob") {
    dt[, delta_prob := data.table::fcoalesce(probability_b, 0) -
                       data.table::fcoalesce(probability_a, 0)]

    return(
      ggplot2::ggplot(dt[!is.na(delta_prob)],
                      ggplot2::aes(x = .plot_class, y = delta_prob,
                                   fill = .plot_class)) +
        ggplot2::geom_violin(alpha = 0.6, trim = TRUE, color = NA) +
        ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.8, alpha = 0.9) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
        ggplot2::scale_fill_manual(values = cols, drop = FALSE) +
        ggplot2::labs(
          title    = "Change in m6A Probability Between Isoforms",
          subtitle = "Undetected sites treated as probability 0",
          x = lab,
          y = expression(Delta * " probability (isoform B - isoform A)"),
          fill = lab) +
        base_theme
    )
  }
}


#' Plot m6A Sites for a Single Gene
#'
#' Track-style plot showing m6A sites on both isoforms of a switch pair.
#'
#' @param gene Character. Gene name or ID.
#' @param m6a_switches data.table from \code{annotate_m6a_switches_genomic()}.
#' @param iso_lengths Optional named numeric vector mapping isoform_id to
#'   length. When supplied, positions are expressed as a fraction [0, 1].
#' @param color_by \code{"isoform_status"} (default) or \code{"m6a_fate"}.
#'
#' @return ggplot object.
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_isoform_details <- function(gene,
                                 m6a_switches,
                                 iso_lengths = NULL,
                                 color_by    = c("isoform_status", "m6a_fate")) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  color_by <- match.arg(color_by)
  if (!color_by %in% names(m6a_switches)) {
    stop("Column '", color_by, "' not found in m6a_switches.")
  }

  .g <- gene
  gene_data <- m6a_switches[gene_id == .g]
  if (nrow(gene_data) == 0) {
    stop(sprintf("No data for gene: %s", gene))
  }

  if (color_by == "isoform_status") {
    lvls <- c("ISOFORM_A_ONLY", "ISOFORM_B_ONLY", "IN_BOTH_ISOFORMS")
    cols <- c(ISOFORM_A_ONLY   = "#D55E00",
              ISOFORM_B_ONLY   = "#009E73",
              IN_BOTH_ISOFORMS = "#0072B2")
    lab  <- "Isoform status"
  } else {
    lvls <- c("LOST", "GAINED", "RETAINED")
    cols <- c(LOST = "#D55E00", GAINED = "#009E73", RETAINED = "#0072B2")
    lab  <- "m6A fate"
  }

  gene_data <- data.table::copy(gene_data)
  gene_data[, .plot_class := factor(get(color_by), levels = lvls)]

  x_col <- if ("start" %in% names(gene_data)) {
    "start"
  } else if ("transcript_position" %in% names(gene_data)) {
    "transcript_position"
  } else {
    "position"
  }

  cond_a_lab <- if (!is.na(gene_data$condition_1[1])) gene_data$condition_1[1] else "condition 1"
  cond_b_lab <- if (!is.na(gene_data$condition_2[1])) gene_data$condition_2[1] else "condition 2"

  dt_a <- gene_data[, .(isoform = isoform_a, position = get(x_col),
                        .plot_class, present = m6a_in_isoform_a,
                        side = paste0("Isoform A (", cond_a_lab, ")"))]
  dt_b <- gene_data[, .(isoform = isoform_b, position = get(x_col),
                        .plot_class, present = m6a_in_isoform_b,
                        side = paste0("Isoform B (", cond_b_lab, ")"))]

  long_dt <- data.table::rbindlist(list(dt_a, dt_b), use.names = TRUE)
  long_dt <- long_dt[present == TRUE]

  if (nrow(long_dt) == 0) {
    stop(sprintf("No m6A-present positions to plot for gene: %s", gene))
  }

  if (!is.null(iso_lengths)) {
    if (is.null(names(iso_lengths))) {
      stop("iso_lengths must be a named vector: names are isoform IDs")
    }
    long_dt[, iso_length := as.numeric(iso_lengths[isoform])]
    long_dt[, position := ifelse(!is.na(iso_length) & iso_length > 0,
                                 position / iso_length, NA_real_)]
    x_lab <- "Relative position within isoform (0 = 5', 1 = 3')"
  } else {
    x_lab <- if (x_col == "start") "Genomic coordinate"
             else if (x_col == "transcript_position") "Transcript position (nt)"
             else "Position"
  }

  subtitle_text <- sprintf("%s \u2192 %s", cond_a_lab, cond_b_lab)
  if (!is.na(gene_data$dif[1]) && !is.na(gene_data$fdr[1])) {
    subtitle_text <- sprintf("dIF = %.2f | FDR = %.2e | %s \u2192 %s",
                             gene_data$dif[1], gene_data$fdr[1],
                             cond_a_lab, cond_b_lab)
  }

  ggplot2::ggplot(long_dt,
                  ggplot2::aes(x = position, y = side, color = .plot_class)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    ggplot2::scale_color_manual(values = cols, drop = FALSE) +
    ggplot2::labs(
      title    = sprintf("m6A Sites in Isoform Switch: %s", gene),
      subtitle = subtitle_text,
      x = x_lab, y = "", color = lab
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "right",
      axis.text.y      = ggplot2::element_text(size = 11, face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold")
    )
}
