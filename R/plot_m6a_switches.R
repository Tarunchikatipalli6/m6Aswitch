#' Plot m6A Switching Events
#'
#' Generates publication-quality plots summarizing m6A changes across isoform switches.
#'
#' @param m6a_switches data.table from annotate_m6a_switches()
#' @param plot_type Type of plot: "summary", "sankey", "heatmap", "upset" (default: "summary")
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
plot_m6a_switches <- function(m6a_switches, plot_type = "summary") {
  
  if (!is.data.table(m6a_switches) || nrow(m6a_switches) == 0) {
    stop("m6a_switches must be a non-empty data.table")
  }
  
  plot_type <- match.arg(plot_type, c("summary", "sankey", "heatmap", "upset"))
  
  if (plot_type == "summary") {
    # Summary: barplot of m6A fates (LOST, GAINED, RETAINED)
    fate_counts <- m6a_switches[, .N, by = m6a_fate]
    
    p <- ggplot2::ggplot(fate_counts, aes(x = m6a_fate, y = N, fill = m6a_fate)) +
      geom_col(position = "dodge") +
      scale_fill_manual(values = c(LOST = "#E74C3C", GAINED = "#27AE60", RETAINED = "#3498DB")) +
      labs(
        title = "m6A Sites Across Isoform Switches",
        x = "m6A Fate",
        y = "Count",
        fill = "Fate"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(p)
    
  } else if (plot_type == "heatmap") {
    # Heatmap: genes x fate
    gene_fate <- m6a_switches[, .N, by = .(gene_id, m6a_fate)]
    
    p <- ggplot2::ggplot(gene_fate, aes(x = m6a_fate, y = gene_id, fill = N)) +
      geom_tile() +
      scale_fill_viridis_c() +
      labs(
        title = "m6A Fate by Gene",
        x = "m6A Fate",
        y = "Gene",
        fill = "Count"
      ) +
      theme_minimal()
    
    return(p)
    
  } else {
    # Placeholder for other plot types
    stop(sprintf("Plot type '%s' not yet implemented.", plot_type))
  }
}

#' Plot Isoform Switch Details
#'
#' Creates a lollipop/track plot showing m6A sites on switching isoforms.
#'
#' @param gene_id Gene to plot
#' @param m6a_switches Annotated m6A switches
#' @param iso_lengths Named vector: isoform_id -> length (for proportional scaling)
#'
#' @return ggplot object
#'
#' @import ggplot2
#' @import data.table
#'
#' @export
plot_isoform_details <- function(gene_id, m6a_switches, iso_lengths = NULL) {
  
  # Filter for gene
  gene_data <- m6a_switches[gene_id == !!gene_id]
  
  if (nrow(gene_data) == 0) {
    stop(sprintf("No data for gene: %s", gene_id))
  }
  
  # Create isoform-level plot
  p <- ggplot2::ggplot(gene_data, aes(x = position, y = factor(isoform_a), color = m6a_fate)) +
    geom_point(size = 3) +
    scale_color_manual(values = c(LOST = "#E74C3C", GAINED = "#27AE60", RETAINED = "#3498DB")) +
    labs(
      title = sprintf("m6A Sites in %s Isoform Switch", gene_id),
      x = "Position",
      y = "Isoform",
      color = "Fate"
    ) +
    theme_minimal() +
    theme(legend.position = "right")
  
  return(p)
}
