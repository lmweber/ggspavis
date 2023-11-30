#' plotFeatureQC
#' 
#' Quality control (QC) plots for spatially resolved transcriptomics data.
#' 
#' Function to generate plots for quality control (QC) purposes for spatially
#' resolved transcriptomics data.
#' 
#' The following types of QC plots are available for per feature QC:
#' 
#' - Histogram (\code{type} = "hist") for a single QC metric, e.g. number of counts
#' per feature. For number of counts per spot, the histogram highlights spots with
#' derived flags, e.g. low abundance genes.
#' - Violin (\code{type} = "violin") for a single QC metric, e.g. number of counts
#' per feature. For number of counts per spot, the violin plot is able to highlights 
#' spots with derived flags, e.g. low abundance genes.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} or \code{SingleCellExperiment} object.
#' 
#' @param type (character) Type of QC plot. Options are "hist", "scatter", and
#'   "spots". See details in description.
#' 
#' @param metric_x (character) Name of column in \code{colData} containing QC
#'   metric to plot on x-axis (e.g. "sum" for number of cells per spot).
#'   Default = "sum". Required for histogram, scatter, and violin plots.
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
#' 
#' @export
#' 
#' @author Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotFeatureQC(spe, type = "hist", metric_x = "cell_count", annotate = "low_abun_gene")
#' plotFeatureQC(spe, type = "scatter", metric_x = "cell_count", metric_y = "sum")

plotFeatureQC <- function(spe, type = c("hist", "violin"), 
                          metric_x = "sum", annotate = NULL, 
                          nbins = 100, pt.size = 0.1,
                          marginal = TRUE, y_reverse = TRUE) {
  
  type <- match.arg(type)
  
  plt_df <- data.frame(rowData(spe))
  
  ## For hist, spot, violin plot
  if (!is.null(annotate)){
    if (!is.logical(plt_df[[annotate]])) {
      stop("For QC, please make sure `annotate` is a binary flag. `metric_x` should be continuous.")
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
  
  
  if (type == "violin") { # must have metric_x (cont), optional annotate (logical)
    
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

