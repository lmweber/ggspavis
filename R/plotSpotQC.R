#' plotSpotQC
#' 
#' Plotting functions for spatial transcriptomics data.
#' 
#' Function to create quality control (QC) plots for spatial transcriptomics
#' data.
#' 
#' The following types of QC plots are available for spot-level or cell-level QC
#' (see \code{\link{plotFeatureQC}} for feature-level QC):
#' 
#' \itemize{
#' \item Histogram (\code{plot_type "histogram"}) for a single QC metric, e.g.
#' number of UMI counts per spot. For number of counts per spot, the histogram
#' can optionally highlight selected spots, e.g. spots with low library size.
#' \item Scatter plot (\code{plot_type = "scatter"}) comparing two QC metrics,
#' e.g. number of detected features vs. number of cells per spot, with optional
#' horizontal and vertical lines highlighting QC filtering thresholds.
#' \item Spot plot (\code{plot_type = "spot"}) showing spots in spatial x-y
#' coordinates, e.g. highlighting selected spots that do not meet filtering
#' thresholds.
#' \item Violin plot (\code{plot_type = "violin"}) for a single QC metric, e.g.
#' number of UMI counts per spot. For number of counts per spot, the violin plot
#' can optionally highlight selected spots, e.g. spots with low library size.
#' }
#' 
#' 
#' @param spe Input data, assumed to be a \code{SpatialExperiment} or
#'   \code{SingleCellExperiment} object.
#' 
#' @param plot_type Type of QC plot. Options are "histogram", "scatter", "spot",
#'   and "violin". See Details for additional details.
#' 
#' @param x_coord Name of column in \code{spatialCoords} (for a
#'   \code{SpatialExperiment} input object) or \code{colData} (for a
#'   \code{SingleCellExperiment} input object) containing x coordinates. Default
#'   = NULL (for a \code{SpatialExperiment}, the first column of
#'   \code{spatialCoords} will be selected in this case). Used for spot plots.
#' 
#' @param y_coord Name of column in \code{spatialCoords} (for a
#'   \code{SpatialExperiment} input object) or \code{colData} (for a
#'   \code{SingleCellExperiment} input object) containing y coordinates. Default
#'   = NULL (for a \code{SpatialExperiment}, the second column of
#'   \code{spatialCoords} will be selected in this case). Used for spot plots.
#' 
#' @param x_metric Name of column in \code{colData} containing QC metric to plot
#'   on x-axis. Required for histograms, scatter plots, and violin plots.
#' 
#' @param y_metric Name of column in \code{colData} containing QC metric to plot
#'   on y-axis. Required for histograms, scatter plots, and violin plots.
#' 
#' @param x_threshold QC filtering threshold on x-axis metric to highlight with
#'   vertical line. Default = NULL. Optional argument used for scatter plots.
#' 
#' @param y_threshold QC filtering threshold on y-axis metric to highlight with
#'   horizontal line. Default = NULL. Optional argument used for scatter plots.
#' 
#' @param trend Whether to show smoothed trend line (loess). Default = TRUE.
#'   Optional argument used for scatter plots.
#' 
#' @param marginal Whether to show marginal histograms. Default = TRUE. Optional
#'   argument used for scatter plots.
#' 
#' @param annotate Name of column in \code{colData} identifying selected spots
#'   that do not meet QC filtering thresholds, which will be highlighted on a
#'   histogram, spot plot, or violin plot. Default = NULL. Optional argument
#'   used for histograms, spot plots, and violin plots.
#' 
#' @param in_tissue Name of column in \code{colData} identifying spots over
#'   tissue (e.g. "in_tissue" for 10x Genomics Visium datasets). If this
#'   argument is provided, only spots over tissue will be shown. Default = NULL.
#'   Optional argument used for spot plots.
#' 
#' @param legend_point_size Legend point size. Default = 3. Optional argument
#'   used for spot plots.
#' 
#' @param n_bins Number of bins for histograms. Default = 100. Optional argument
#'   used for histograms.
#' 
#' @param point_size Point size. Default = 0.3. Optional argument for scatter
#'   plots, spot plots, and violin plots. Suggested values: 0.5 for scatter
#'   plots, 0.3 for spot plots, 0.1 for violin plots.
#' 
#' @param y_reverse Whether to reverse y coordinates. This is usually required
#'   for 10x Genomics Visium datasets when using the default coordinate values.
#'   Default = TRUE. Set to FALSE if not needed, e.g. for other platforms.
#'   Optional argument used for spot plots.
#' 
#' 
#' @return Returns a ggplot object, which may be further modified using ggplot
#'   functions.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom ggplot2 ggplot aes_string geom_histogram geom_point geom_vline
#'   geom_hline geom_smooth geom_violin geom_jitter scale_fill_manual
#'   scale_color_manual xlab ylab labs coord_fixed theme_bw theme
#'   element_text element_blank guides scale_y_reverse
#' @importFrom ggside geom_xsidehistogram geom_ysidehistogram
#' 
#' 
#' @export
#' 
#' @author Lukas M. Weber and Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' 
#' spe$sum <- colSums(counts(spe))
#' spe$low_libsize <- spe$sum < 400
#' 
#' plotSpotQC(spe, plot_type = "histogram", x_metric = "sum", annotate = "low_libsize")
#' plotSpotQC(spe, plot_type = "scatter", x_metric = "sum", y_metric = "cell_count")
#' plotSpotQC(spe, plot_type = "spot", annotate = "low_libsize")
#' plotSpotQC(spe, plot_type = "violin", x_metric = "sum", annotate = "low_libsize")
#' 
plotSpotQC <- function(spe, 
                       plot_type = c("histogram", "scatter", "spot", "violin"), 
                       x_coord = NULL, y_coord = NULL, 
                       x_metric = NULL, y_metric = NULL, 
                       x_threshold = NULL, y_threshold = NULL, 
                       trend = TRUE, marginal = TRUE, 
                       annotate = NULL, in_tissue = NULL, 
                       legend_point_size = 3, 
                       n_bins = 100, point_size = 0.3, 
                       y_reverse = TRUE) {
  
  # check validity of arguments
  plot_type <- match.arg(plot_type)
  
  if (!is.null(in_tissue)) {
    stopifnot(is.character(in_tissue))
  }
  
  if (!is.null(annotate)) {
    stopifnot(is.character(annotate))
    if (!(annotate %in% colnames(colData(spe)))) {
      stop("'annotate' should be the name of a column in colData")
    }
  }
  
  # set up data frame for plotting
  if (class(spe) == "SpatialExperiment") {
    # select default columns of x and y coordinates
    if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
    if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
    df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  } else if (class(spe) == "SingleCellExperiment") {
    if (is.null(x_coord) || is.null(y_coord)) {
      stop("Please provide 'x_coord' and 'y_coord' arguments to specify ", 
           "columns in colData containing x and y coordinates.")
    }
    df <- as.data.frame(colData(spe))
  }
  
  # for histogram, spot, or violin plots
  if (!is.null(annotate)) {
    if (!is.logical(df[[annotate]])) {
      stop("For QC plots, 'annotate' should be a logical vector, and 'x_metric' ", 
           "and/or 'y_metric' should be vectors of continuous values.")
    }
  }
  
  if (!is.null(in_tissue)) {
    df <- df[df[, in_tissue] == 1, ]
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
    
    p <- p + theme_bw()
  }
  
  
  # scatter: requires 'x_metric' (continuous) and 'y_metric' (continuous),
  # additional optional arguments
  
  if (plot_type == "scatter") {
    
    p <- ggplot(df, aes_string(x = x_metric, y = y_metric)) + 
      geom_point(size = point_size) + 
      theme_bw()
    
    if (!is.null(x_threshold)) {
      p <- p + 
        geom_vline(xintercept = x_threshold, color = "red")
    }
    if (!is.null(y_threshold)) {
      p <- p + 
        geom_hline(yintercept = y_threshold, color = "red")
    }
    if (trend) {
      p <- p + 
        geom_smooth(method = "loess", se = FALSE)
    }
    if (marginal) {
      p <- p + 
        geom_xsidehistogram() + 
        geom_ysidehistogram()
    }
  }
  
  
  # spot: requires 'x_coord' (continuous), 'y_coord' (continuous), optionally
  # 'annotate' (logical)
  
  if (plot_type == "spot") {
    
    # spots at 'x_coord' and 'y_coord', optionally colored by 'annotate'
    if (!is.null(annotate)) {
      p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = annotate)) + 
        geom_point(size = point_size) + 
        scale_color_manual(values = c("gray85", "red")) + 
        guides(color = guide_legend(override.aes = list(size = legend_point_size)))
    } else if (is.null(annotate)) {
      p <- ggplot(df, aes_string(x = x_coord, y = y_coord)) + 
        geom_point(size = point_size) + 
        scale_color_manual(values = "gray85")
    }
    
    p <- p + 
      coord_fixed() + 
      theme_bw() + 
      theme(panel.grid = element_blank(), 
            axis.title = element_blank(), 
            axis.text = element_blank(), 
            axis.ticks = element_blank())
    
    if (y_reverse) {
      p <- p + scale_y_reverse()
    }
  }
  
  
  # violin: requires 'x_metric' (continuous), optionally 'annotate' (logical)
  
  if (plot_type == "violin") {
    
    df[["dummy"]] <- rep(" ", nrow(df))
    
    p <- ggplot(df, aes_string(x = "dummy", y = x_metric, fill = "dummy")) + 
      geom_violin(trim = TRUE, alpha = 0.9) + 
      scale_fill_manual(values = c("gray70")) + 
      xlab("Sample") + 
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
  }
  
  # return plot
  p
}
