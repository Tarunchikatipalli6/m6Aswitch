#' Plot m6A Switching Events
#'
#' Generates publication-quality plots summarizing m6A changes across isoform switches.
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param plot_type Type of plot:
#'   "summary", "heatmap", "delta_prob", "motif_enrichment", "metagene"
#'   (default: "summary")
#' @param top_n_genes Number of top genes to display in heatmap (default: 30)
#'
#' @return ggplot object
#'
#' @examples
#' \dontrun{
#' plot_m6a_switches(m6a_switches, plot_type = "summary")
#' }
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_m6a_switches <- function(m6a_switches,
                              plot_type = "summary",
                              top_n_genes = 30) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  plot_type <- match.arg(plot_type, c(
    "summary", "heatmap", "delta_prob", "motif_enrichment", "metagene"
  ))

  fate_levels <- c("LOST", "GAINED", "RETAINED")
  fate_colors <- c(
    LOST = "#D55E00",
    GAINED = "#009E73",
    RETAINED = "#0072B2"
  )

  dt <- data.table::copy(m6a_switches)
  dt[, m6a_fate := factor(m6a_fate, levels = fate_levels)]

  base_theme <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  if (plot_type == "summary") {
    fate_counts <- dt[, .N, by = m6a_fate][order(m6a_fate)]
    fate_counts[, pct := 100 * N / sum(N)]

    p <- ggplot2::ggplot(fate_counts, ggplot2::aes(x = m6a_fate, y = N, fill = m6a_fate)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%d\n(%.1f%%)", N, pct)),
        vjust = -0.4,
        size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
      ggplot2::labs(
        title = "m6A Site Fate Across Isoform Switches",
        x = "m6A fate",
        y = "Number of sites",
        fill = "Fate"
      ) +
      base_theme

    return(p)
  }

  if (plot_type == "heatmap") {
    gene_totals <- dt[, .N, by = gene_id][order(-N)]
    keep_genes <- head(gene_totals$gene_id, top_n_genes)

    gene_fate <- dt[gene_id %in% keep_genes, .N, by = .(gene_id, m6a_fate)]
    gene_fate[, gene_id := factor(gene_id, levels = rev(keep_genes))]

    p <- ggplot2::ggplot(gene_fate, ggplot2::aes(x = m6a_fate, y = gene_id, fill = N)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.2) +
      ggplot2::scale_fill_viridis_c(option = "C", name = "Count") +
      ggplot2::labs(
        title = sprintf("m6A Fate by Gene (Top %d genes)", length(keep_genes)),
        x = "m6A fate",
        y = "Gene"
      ) +
      base_theme

    return(p)
  }

  if (plot_type == "delta_prob") {
    # numeric shift in methylation strength where possible
    dt[, delta_prob := probability_b - probability_a]

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
        x = "m6A fate",
        y = expression(Delta * " probability (isoform B - isoform A)"),
        fill = "Fate"
      ) +
      base_theme

    return(p)
  }

  if (plot_type == "motif_enrichment") {
    if (!all(c("drach_motif_a", "drach_motif_b") %in% names(dt))) {
      stop("motif_enrichment requires columns: drach_motif_a and drach_motif_b")
    }

    dt[, drach_any := as.logical(drach_motif_a) | as.logical(drach_motif_b)]
    motif_summary <- dt[, .(
      total = .N,
      drach_true = sum(drach_any, na.rm = TRUE)
    ), by = m6a_fate]
    motif_summary[, prop := drach_true / total]

    p <- ggplot2::ggplot(motif_summary, ggplot2::aes(x = m6a_fate, y = prop, fill = m6a_fate)) +
      ggplot2::geom_col(width = 0.75) +
      ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.1f%%", 100 * prop)),
        vjust = -0.4,
        size = 3.5
      ) +
      ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
      ggplot2::scale_y_continuous(limits = c(0, 1), labels = function(x) paste0(round(100 * x), "%")) +
      ggplot2::labs(
        title = "DRACH Motif Enrichment by m6A Fate",
        x = "m6A fate",
        y = "Fraction DRACH-positive",
        fill = "Fate"
      ) +
      base_theme

    return(p)
  }

  if (plot_type == "metagene") {
    if (!"relative_position" %in% names(dt)) {
      stop("metagene plot requires column 'relative_position' in [0,1].")
    }

    p <- ggplot2::ggplot(
      dt[!is.na(relative_position)],
      ggplot2::aes(x = relative_position, color = m6a_fate, fill = m6a_fate)
    ) +
      ggplot2::geom_density(alpha = 0.2, linewidth = 1) +
      ggplot2::scale_color_manual(values = fate_colors, drop = FALSE) +
      ggplot2::scale_fill_manual(values = fate_colors, drop = FALSE) +
      ggplot2::labs(
        title = "Metagene Distribution of m6A Sites",
        x = "Relative transcript position (0 = 5', 1 = 3')",
        y = "Density",
        color = "Fate",
        fill = "Fate"
      ) +
      base_theme

    return(p)
  }

  stop(sprintf("Unsupported plot_type: %s", plot_type))
}

#' Plot Isoform Switch Details
#'
#' Creates a lollipop/track-style plot showing m6A sites on both switching isoforms.
#'
#' @param gene_id Gene to plot
#' @param m6a_switches Annotated m6A switches
#' @param iso_lengths Optional named numeric vector: isoform_id -> length
#'
#' @return ggplot object
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_isoform_details <- function(gene_id, m6a_switches, iso_lengths = NULL) {

  if (!data.table::is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }

  # FIX: avoid data.table NSE conflict with function argument name
  gene_data <- m6a_switches[get("gene_id") == gene_id]

  if (nrow(gene_data) == 0) {
    stop(sprintf("No data for gene: %s", gene_id))
  }

  fate_colors <- c(
    LOST = "#D55E00",
    GAINED = "#009E73",
    RETAINED = "#0072B2"
  )

  # Build long table so both isoform_a and isoform_b are plotted
  dt_a <- gene_data[, .(
    isoform = isoform_a,
    position = position,
    m6a_fate = m6a_fate,
    present = m6a_in_isoform_a,
    side = "A"
  )]
  dt_b <- gene_data[, .(
    isoform = isoform_b,
    position = position,
    m6a_fate = m6a_fate,
    present = m6a_in_isoform_b,
    side = "B"
  )]

  long_dt <- data.table::rbindlist(list(dt_a, dt_b), use.names = TRUE)
  long_dt <- long_dt[present == TRUE]

  if (nrow(long_dt) == 0) {
    stop(sprintf("No m6A-present positions to plot for gene: %s", gene_id))
  }

  # optional scaling
  if (!is.null(iso_lengths)) {
    if (is.null(names(iso_lengths))) stop("iso_lengths must be a named vector: names are isoform IDs")
    long_dt[, iso_length := as.numeric(iso_lengths[isoform])]
    long_dt[, position_scaled := ifelse(!is.na(iso_length) & iso_length > 0, position / iso_length, NA_real_)]
    x_var <- "position_scaled"
    x_lab <- "Relative position within isoform"
  } else {
    x_var <- "position"
    x_lab <- "Position"
  }

  p <- ggplot2::ggplot(
    long_dt,
    ggplot2::aes(x = .data[[x_var]], y = factor(isoform), color = m6a_fate, shape = side)
  ) +
    ggplot2::geom_point(size = 3, alpha = 0.9) +
    ggplot2::scale_color_manual(values = fate_colors, drop = FALSE) +
    ggplot2::labs(
      title = sprintf("m6A Sites in Isoform Switches: %s", gene_id),
      x = x_lab,
      y = "Isoform",
      color = "Fate",
      shape = "Switch side"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "right"
    )

  return(p)
}
