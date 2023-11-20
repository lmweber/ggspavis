#' plotQC
#' 
#' Quality control (QC) plots for spatially resolved transcriptomics data.
#' 
#' Function to generate plots for quality control (QC) purposes for spatially
#' resolved transcriptomics data.
#' 
#' The following types of QC plots are available:
#' 
#' - Histogram (\code{type} = "bar") for a single QC metric, e.g. number of counts
#' per spot. For number of counts per spot, the histogram highlights spots with
#' derived flags, e.g. low library size.
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
#' @param type (character) Type of QC plot. Options are "hist", "scatter", and
#'   "spots". See details in description.
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates. Default = NULL, which selects the first column. Used for
#'   spot-based plots.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates. Default = NULL, which selects the second column. Used for
#'   spot-based plots.
#' 
#' @param in_tissue (character) Name of column in \code{colData} identifying
#'   spots over tissue, e.g. "in_tissue" for 10x Genomics Visium data. If this
#'   argument is provided, only spots over tissue will be shown. Alternatively,
#'   set to NULL to display all spots. Default = "in_tissue".
#' 
#' @param metric_x (character) Name of column in \code{colData} containing QC
#'   metric to plot on x-axis (e.g. "cell_count" for number of cells per spot).
#'   Default = "cell_count". Required for histogram and scatterplots.
#' 
#' @param metric_y (character) Name of column in \code{colData} containing QC
#'   metric to plot on y-axis (e.g. "sum" for number of detected transcripts, or
#'   "high_mito" for binary ). Default = "sum". Required for
#'   histogram and scatterplots.
#' 
#' @param discard (character) Name of column in \code{colData} identifying
#'   discarded spots that do not meet filtering thresholds, which will be
#'   highlighted on a spot-based plot. Default = "discard". Optional for
#'   spot-based plots.
#'   
#' @param nbins (numeric) Adjusting the histogram bin width. Optional for
#'   spot-based and scatterplots.
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
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes_string geom_histogram geom_point geom_smooth 
#'   geom_hline geom_vline coord_fixed labs ggtitle theme_bw theme element_blank 
#'   scale_y_reverse scale_fill_manual scale_color_manual
#' @importFrom ggside geom_xsidehistogram geom_ysidehistogram
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotQC(spe, type = "hist", metric_x = "cell_count")
#' colData(spe)$sum <- colSums(counts(spe))
#' plotQC(spe, type = "scatter", metric_x = "cell_count", metric_y = "sum")
#' 
plotQC <- function(spe, type = c("hist", "scatter", "spots"), 
                   x_coord = NULL, y_coord = NULL, in_tissue = "in_tissue", 
                   metric_x = "cell_count", metric_y = "sum", 
                   discard = "discard", nbins = 100, 
                   threshold_x = NULL, threshold_y = NULL, 
                   trend = TRUE, marginal = TRUE, y_reverse = TRUE) {
  
  type <- match.arg(type)
  if (!is.null(in_tissue)) stopifnot(is.character(in_tissue))
  
  if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
  
  if(class(spe) == "SingleCellExperiment"){
    plt_df <- cbind.data.frame(colData(spe))
  }else if(class(spe) == "SpatialExperiment"){
    plt_df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  }
  
  if (type == "hist") {
    stopifnot(is.logical(plt_df[, metric_y]))
    stopifnot(is.numeric(nbins))
    
    p <- ggplot(plt_df, aes_string(x = metric_x, fill = metric_y)) +
      geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = nbins) +
      scale_fill_manual(values = c("gray70", "firebrick2")) +
      theme_bw() +
      labs(fill = metric_y) + xlab(metric_x) +
      theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12)) + 
      ggtitle("QC metrics")
  }
  
  if (type == "scatter") {
    p <- ggplot(plt_df, aes_string(x = metric_x, y = metric_y)) + 
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
    
    if (!is.null(in_tissue)) {
      plt_df <- plt_df[plt_df[, in_tissue] == 1, ]
    }
    
    p <- ggplot(plt_df, aes_string(x = x_coord, y = y_coord, color = discard)) + 
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

