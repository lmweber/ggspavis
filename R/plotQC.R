#' plotQC
#' 
#' Quality control (QC) plots for spatially resolved transcriptomics data.
#' 
#' Function to generate plots for quality control (QC) purposes for spatially
#' resolved transcriptomics data.
#' 
#' The following types of QC plots are available:
#' 
#' - Barplot (\code{type} = "bar") for a single QC metric, e.g. number of cells
#' per spot. For number of cells per spot, the barplot highlights spots with
#' zero cells.
#' - Scatterplot (\code{type} = "scatter") comparing two QC metrics, e.g. number
#' of detected features vs. number of cells per spot, with optional vertical and
#' horizontal lines highlighting QC filtering thresholds.
#' - Spots (\code{type} = "spots") i.e. spots in spatial (x-y) coordinates,
#' highlighting discarded spots that do not meet filtering thresholds.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param type (character) Type of QC plot. Options are "bar", "scatter", and
#'   "spots". See details in description.
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates. Default = "x". Required for spot-based plots.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates. Default = "y". Required for spot-based plots.
#' 
#' @param in_tissue (logical) Whether to show only spots over tissue, or all
#'   spots. Options are TRUE (show spots over tissue; requires a column labelled
#'   "in_tissue" in \code{spatialData} identifying spots over tissue, which is
#'   the standard format for 10x Genomics Visium data), FALSE (show all spots),
#'   or a character value with the name of a column in \code{spatialData}
#'   identifying the spots to show.
#' 
#' @param metric_x (character) Name of column in \code{colData} containing QC
#'   metric to plot on x-axis (e.g. "cell_count" for number of cells per spot).
#'   Default = "cell_count". Required for barplots and scatterplots.
#' 
#' @param metric_y (character) Name of column in \code{colData} containing QC
#'   metric to plot on y-axis (e.g. "sum" for number of detected transcripts, or
#'   "detected" for number of detected genes). Default = "sum". Required for
#'   scatterplots.
#' 
#' @param discard (character) Name of column in \code{colData} identifying
#'   discarded spots that do not meet filtering thresholds, which will be
#'   highlighted on a spot-based plot. Default = "discard". Optional for
#'   spot-based plots.
#' 
#' @param highlight_zeros (logical) Whether to highlight bar for x = 0 (e.g.
#'   zero cells per spot). Default = TRUE. Optional for barplots.
#' 
#' @param threshold_x (numeric) Filtering threshold for QC metric on x-axis,
#'   which will be highlighted with a vertical bar. Default = NULL. Optional for
#'   scatterplots.
#' 
#' @param threshold_y (numeric) Filtering threshold for QC metric on y-axis,
#'   which will be highlighted with a horizontal bar. Default = NULL. Optional
#'   for scatterplots.
#' 
#' @param trend (logical) Whether to display a smoothed trend (loess) for
#'   scatterplots. Default = TRUE. Optional for scatterplots.
#' 
#' @param marginal (logical) Whether to display marginal histograms for
#'   scatterplots. Default = TRUE. Optional for scatterplots.
#' 
#' @param y_reverse (logical) Whether to reverse y coordinates, which is often
#'   required for 10x Genomics Visium data. Default = TRUE.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialData spatialCoords
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes_string geom_bar geom_point geom_smooth 
#'   geom_hline geom_vline coord_fixed labs ggtitle theme_bw theme element_blank 
#'   scale_y_reverse scale_fill_manual scale_color_manual
#' @importFrom ggside geom_xsidehistogram geom_ysidehistogram
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotQC(spe, type = "bar", metric_x = "cell_count")
#' colData(spe)$sum <- colSums(counts(spe))
#' plotQC(spe, type = "scatter", metric_x = "cell_count", metric_y = "sum")
#' 
plotQC <- function(spe, type = c("bar", "scatter", "spots"), 
                   x_coord = "x", y_coord = "y", in_tissue = TRUE, 
                   metric_x = "cell_count", metric_y = "sum", 
                   discard = "discard", highlight_zeros = TRUE, 
                   threshold_x = NULL, threshold_y = NULL, 
                   trend = TRUE, marginal = TRUE, y_reverse = TRUE) {
  
  type <- match.arg(type)
  
  stopifnot(is.character(x_coord) & is.character(y_coord))
  
  df <- as.data.frame(cbind(colData(spe), spatialData(spe), spatialCoords(spe)))
  
  if (type == "bar") {
    p <- ggplot(df, aes_string(x = metric_x)) + 
      geom_bar(fill = "gray70") + 
      labs(x = metric_x) + 
      ggtitle("QC metrics") + 
      theme_bw()
    
    if (highlight_zeros) {
      df$is_zero <- factor(df[, metric_x] == 0)
      p <- p + 
        geom_bar(data = df, aes_string(fill = "is_zero")) + 
        scale_fill_manual(values = c("gray70", "firebrick2"))
    }
  }
  
  if (type == "scatter") {
    p <- ggplot(df, aes_string(x = metric_x, y = metric_y)) + 
      geom_point(size = 0.5) + 
      ggtitle("QC metrics") + 
      theme_bw()
    
    if (!is.null(threshold_x)) {
      p <- p + geom_vline(xintercept = threshold_x, color = "red")
    }
    if (!is.null(threshold_y)) {
      p <- p + geom_hline(yintercept = threshold_y, color = "red")
    }
    if (trend) {
      p <- p + geom_smooth(method = "loess", se = FALSE)
    }
    if (marginal) {
      p <- p + 
        geom_xsidehistogram() + 
        geom_ysidehistogram()
    }
  }
  
  if (type == "spots") {
    p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = discard)) + 
      geom_point(size = 0.3) + 
      coord_fixed() + 
      scale_color_manual(values = c("gray85", "red")) + 
      ggtitle("QC spots") + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (y_reverse) {
      p <- p + scale_y_reverse()
    }
  }
  
  p
}

