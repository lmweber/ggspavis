#' plotFeatureQC
#' 
#' Plotting functions for spatial transcriptomics data.
#' 
#' Function to create quality control (QC) plots for spatial transcriptomics
#' data.
#' 
#' The following types of QC plots are available for feature-level QC (see
#' \code{\link{plotSpotQC}} for spot-level or cell-level QC):
#' 
#' \itemize{
#' \item Histogram (\code{plot_type = "histogram"}) for a single QC metric, e.g.
#' total UMI counts across all spots per feature. The histogram can optionally
#' highlight selected features, e.g. low abundance features.
#' \item Violin (\code{plot_type = "violin"}) for a single QC metric, e.g. total
#' UMI counts across all spots per feature. The violin plot can optionally
#' highlight selected features, e.g. low abundance features.
#' }
#' 
#' 
#' @param spe Input data, assumed to be a \code{SpatialExperiment} or
#'   \code{SingleCellExperiment} object.
#' 
#' @param plot_type Type of QC plot. Options are "histogram" and "violin". See
#'   Details for additional details.
#' 
#' @param x_metric Name of column in \code{rowData} containing feature-level QC
#'   metric to plot on x-axis. Required for histograms and violin plots.
#' 
#' @param annotate Name of column in \code{rowData} identifying selected
#'   features that do not meet QC filtering thresholds, which will be
#'   highlighted on a histogram or violin plot. Default = NULL. Optional
#'   argument used for histograms and violin plots.
#' 
#' @param n_bins Number of bins for histograms. Default = 100. Optional argument
#'   used for histograms.
#' 
#' @param point_size Point size. Default = 0.1. Optional argument for violin
#'   plots.
#' 
#' @param scale_log1p Whether to log1p-scale axes. Default = TRUE.
#' 
#' 
#' @return Returns a ggplot object, which may be further modified using ggplot
#'   functions.
#' 
#' 
#' @importFrom SummarizedExperiment rowData
#' @importFrom ggplot2 ggplot aes_string geom_histogram geom_violin geom_jitter
#'   scale_fill_manual scale_color_manual xlab labs theme_bw theme element_blank
#' 
#' 
#' @export
#' 
#' @author Yixing E. Dong and Lukas M. Weber
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' 
#' rowData(spe)$feature_sum <- rowSums(counts(spe))
#' rowData(spe)$low_abundance <- rowData(spe)$feature_sum < 100
#' 
#' plotFeatureQC(spe, plot_type = "histogram", 
#'               x_metric = "feature_sum", annotate = "low_abundance")
#' plotFeatureQC(spe, plot_type = "violin", 
#'               x_metric = "feature_sum", annotate = "low_abundance")
#' 
plotFeatureQC <- function(spe, plot_type = c("histogram", "violin"), 
                          x_metric = NULL, annotate = NULL, 
                          n_bins = 100, point_size = 0.1, 
                          scale_log1p = TRUE) {
  
  # check validity of arguments
  plot_type <- match.arg(plot_type)
  
  df <- data.frame(rowData(spe))
  
  # for histogram or violin plots
  if (!is.null(annotate)) {
    stopifnot(is.character(annotate))
    if (!is.logical(df[[annotate]])) {
      stop("'annotate' should be a vector of logical values in rowData")
    }
  }
  
  
  # histogram: requires 'x_metric' (continuous), optionally 'annotate' (logical)
  
  if (plot_type == "histogram") {
    
    stopifnot(is.numeric(n_bins))
    
    # histogram showing 'x_metric', optionally colored by 'annotate'
    if (!is.null(annotate)) {
      p <- ggplot(df, aes_string(x = x_metric, fill = annotate)) + 
        geom_histogram(bins = n_bins, color = "#e9ecef", alpha = 0.6, 
                       position = "identity") + 
        scale_fill_manual(values = c("gray70", "red")) + 
        xlab(x_metric) + 
        labs(fill = annotate)
    } else if (is.null(annotate)) {
      p <- ggplot(df, aes_string(x = x_metric)) + 
        geom_histogram(bins = n_bins, color = "#e9ecef", alpha = 0.6, 
                       position = "identity") + 
        scale_fill_manual(values = c("gray70")) + 
        xlab(x_metric)
    }
    
    if (scale_log1p) {
      p <- p + scale_x_continuous(transform = "log1p")
    }
    
    p <- p + theme_bw()
  }
  
  
  # violin: requires 'x_metric' (continuous), optionally 'annotate' (logical)
  
  if (plot_type == "violin") {
    
    df[["dummy"]] <- rep(" ", nrow(df))
    
    p <- ggplot(df, aes_string(x = "dummy", y = x_metric, fill = "dummy")) + 
      geom_violin(trim = TRUE, alpha = 0.9) + 
      scale_fill_manual(values = c("gray70")) + 
      ylab(x_metric) + 
      theme_bw() + 
      theme(legend.position="none", 
            panel.grid = element_blank())
    
    if (is.null(annotate)) {
      # violins for 'x_metric'
      p <- p + 
        geom_jitter(size = point_size)
    } else if (!is.null(annotate)) {
      # violins for 'x_metric', colored by 'annotate' (colors in order FALSE, TRUE)
      p <- p + 
        geom_jitter(aes_string(color = annotate), size = point_size) + 
        scale_color_manual(values = c("black", "red"))
    }
    
    if (scale_log1p) {
      p <- p + scale_y_continuous(transform = "log1p")
    }
  }
  
  # return plot
  p
}
