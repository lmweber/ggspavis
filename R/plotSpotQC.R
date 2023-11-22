#' plotSpotQC
#' 
#' Quality control (QC) plots for spatially resolved transcriptomics data.
#' 
#' Function to generate plots for quality control (QC) purposes for spatially
#' resolved transcriptomics data.
#' 
#' The following types of QC plots are available for per spot/cell QC:
#' 
#' - Histogram (\code{type} = "hist") for a single QC metric, e.g. number of counts
#' per spot. For number of counts per spot, the histogram highlights spots with
#' derived flags, e.g. low library size.
#' - Scatterplot (\code{type} = "scatter") comparing two QC metrics, e.g. number
#' of detected features vs. number of cells per spot, with optional vertical and
#' horizontal lines highlighting QC filtering thresholds.
#' - Spots (\code{type} = "spots") i.e. spots in spatial (x-y) coordinates,
#' highlighting flagged spots that do not meet filtering thresholds.
#' - Violin (\code{type} = "violin") for a single QC metric, e.g. number of counts
#' per spot. For number of counts per spot, the violin plot is able to highlights 
#' spots with derived flags, e.g. low library size.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} or \code{SingleCellExperiment} object.
#' 
#' @param type (character) Type of QC plot. Options are "hist", "scatter",
#'   "spots", and "violin". See details in description.
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
#'   set to NULL to display all spots for Xenium, CosMx, etc. Default = "in_tissue".
#' 
#' @param metric_x (character) Name of column in \code{colData} containing QC
#'   metric to plot on x-axis (e.g. "sum" for number of cells per spot).
#'   Default = "sum". Required for histogram, scatter, and violin plots.
#' 
#' @param metric_y (character) Name of column in \code{colData} containing QC
#'   metric to plot on y-axis (e.g. "cell_count" for number of detected transcripts, or
#'   "high_mito" for binary ). Default = "cell_count". Only required for scatter plots.
#' 
#' @param annotate (logical) Name of column in \code{colData} identifying
#'   flagged spots that do not meet filtering thresholds, which will be
#'   highlighted on a histogram, spot-based, or violin plot. Default = NULL. 
#'   Not needed for scatter plot.
#'   
#' @param nbins (numeric) Adjusting the histogram bin width. Optional for
#'   spot-based and scatter and violin plots.
#'   
#' @param pt.size (numeric) Adjusting the scatter, spot-based, violin jitter point size. 
#'   Optional for histogram. Suggested point size for scatter = 0.5, for spot-based
#'   = 0.3, and for violin jitter = 0.1.
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
#' spe$sum <- colSums(counts(spe))
#' spe$low_libsize <- spe$sum < 400
#' 
#' plotSpotQC(spe, type = "hist", metric_x = "sum", annotate = "low_libsize")
#' plotSpotQC(spe, type = "scatter", metric_x = "sum", metric_y = "cell_count")
#' plotSpotQC(spe, type = "spots", annotate = "low_libsize")
#' plotSpotQC(spe, type = "violin", metric_x = "sum", annotate = "low_libsize")

plotSpotQC <- function(spe, type = c("hist", "scatter", "spots", "violin"), 
                   x_coord = NULL, y_coord = NULL, in_tissue = "in_tissue", 
                   metric_x = "sum", metric_y = "cell_count", 
                   annotate = NULL, nbins = 100, pt.size = 0.3,
                   threshold_x = NULL, threshold_y = NULL, 
                   trend = TRUE, marginal = TRUE, y_reverse = TRUE) {
  
  type <- match.arg(type)
  if (!is.null(in_tissue)) stopifnot(is.character(in_tissue))
  
  if(class(spe) == "SingleCellExperiment"){
    if (is.null(x_coord)) {stop("Please specify x_coord name in the SCE colData.")}
    if (is.null(y_coord)) {stop("Please specify y_coord name in the SCE colData.")}
    
    plt_df <- cbind.data.frame(colData(spe))
  }else if(class(spe) == "SpatialExperiment"){
    if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
    if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
    
    plt_df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  }
  
  ## For hist, spot, violin plot
  if (!is.null(annotate)){
    if (!is.logical(plt_df[[annotate]])) {
    stop("For QC, please make sure `annotate` is a binary flag. `metric_x` or/and `metric_y` should be continuous.")
    }
  }
  
  
  if (type == "hist") { # must have metric_x (cont), optional annotate (logical)
    stopifnot(is.numeric(nbins))
    
    if (!is.null(annotate)){ # Histogram of metric_x (cont) , colored by annotate (logical)
      p <- ggplot(plt_df, aes_string(x = metric_x, fill = annotate)) +
        geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = nbins) +
        scale_fill_manual(values = c("gray70", "red")) +
        labs(fill = annotate) + xlab(metric_x)
    }else if(is.null(annotate)){ # Just gray histogram of metric_x (cont) 
      p <- ggplot(plt_df, aes_string(x = metric_x)) +
        geom_histogram(color = "#e9ecef", alpha = 0.6, position = 'identity', bins = nbins) +
        scale_fill_manual(values = c("gray70")) +
        xlab(metric_x)
    }
    
    p <- p +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "plain", size = 12)) + 
      ggtitle("QC metrics")
  }
  
  if (type == "scatter") { # must have metric_x (cont), metric_x (cont)
    p <- ggplot(plt_df, aes_string(x = metric_x, y = metric_y)) + 
      geom_point(size = pt.size) + 
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
  
  if (type == "spots") { # must have x_coord (cont), y_coord (cont), optional annotate (logical)
    
    if (!is.null(in_tissue)) {
      plt_df <- plt_df[plt_df[, in_tissue] == 1, ]
    }
    
    if (!is.null(annotate)){ # At x_coord (cont) & y_coord (cont), colored by annotate (logical)
      p <- ggplot(plt_df, aes_string(x = x_coord, y = y_coord, color = annotate)) + 
        geom_point(size = pt.size) + 
        scale_color_manual(values = c("gray85", "red"))
    }else if(is.null(annotate)){ # just gray spots at x_coord (cont) & y_coord (cont)
      p <- ggplot(plt_df, aes_string(x = x_coord, y = y_coord)) + 
        geom_point(size = pt.size) + 
        scale_color_manual(values = "gray85")
    }
    
    p <- p + 
      coord_fixed() + 
      ggtitle("QC spots") + 
      theme_bw() + 
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())
    
    if (y_reverse) {
      p <- p + scale_y_reverse()
    }
  }
  
  if (type == "violin") { # must have metric_x (cont), optional annotate (logical)
    
    if (!is.null(in_tissue)) {
      plt_df <- plt_df[plt_df[, in_tissue] == 1, ]
    }
    
    plt_df$dummy <- rep(" ", nrow(plt_df))
  
    p <- ggplot(plt_df, aes_string(x="dummy", y=metric_x, fill="dummy")) +
      geom_violin(trim=TRUE, alpha = 0.9) + 
      scale_fill_manual(values = c("gray70")) +
      theme_bw() + 
      xlab("Sample") + ylab("") + ggtitle(metric_x) + 
      theme(legend.position="none", 
            panel.grid = element_blank(), 
            panel.border = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
            axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")
            ) 
    
    if(is.null(annotate)){
      p <- p + geom_jitter(size = pt.size) # just violin for metric_x (cont)
    }else if(!is.null(annotate)){
      p <- p + geom_jitter(aes_string(color = annotate), size = pt.size) + # violin for metric_x (cont), colored by annotate (logical)
        scale_color_manual(values = c("black", "red")) # Note: in order of FALSE, TRUE
    }
    
     
  }
  
  p
}

