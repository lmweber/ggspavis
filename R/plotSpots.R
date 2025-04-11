#' plotSpots
#' 
#' Plotting functions for spatial transcriptomics data.
#' 
#' Function to create spot plot showing spatial locations in x-y coordinates
#' with optional annotations such as cluster labels, expression of a gene, or
#' quality control metrics.
#' 
#' 
#' @param spe Input data, assumed to be a \code{SpatialExperiment} or
#'   \code{SingleCellExperiment} object.
#' 
#' @param x_coord Name of column in \code{spatialCoords} (for a
#'   \code{SpatialExperiment} input object) or \code{colData} (for a
#'   \code{SingleCellExperiment} input object) containing x coordinates. Default
#'   = NULL (for a \code{SpatialExperiment}, the first column of
#'   \code{spatialCoords} will be selected in this case).
#' 
#' @param y_coord Name of column in \code{spatialCoords} (for a
#'   \code{SpatialExperiment} input object) or \code{colData} (for a
#'   \code{SingleCellExperiment} input object) containing y coordinates. Default
#'   = NULL (for a \code{SpatialExperiment}, the second column of
#'   \code{spatialCoords} will be selected in this case).
#' 
#' @param sample_id Name of column in \code{colData} containing sample IDs. This
#'   argument is only required for datasets containing multiple samples (tissue
#'   sections). If provided, samples will be shown in multiple panels using
#'   facetting. Default = NULL.
#' 
#' @param in_tissue Name of column in \code{colData} identifying spots over
#'   tissue (e.g. "in_tissue" for 10x Genomics Visium datasets). If this
#'   argument is provided, only spots over tissue will be shown. Default =
#'   "in_tissue". Set to NULL to display all spots.
#' 
#' @param annotate Variable to show as annotations. This may be discrete or
#'   continuous. For a discrete variable (e.g. cluster labels), this should be
#'   the name of a column in \code{colData} containing a character vector or
#'   factor. For a continuous variable (e.g. a gene name), this should be an
#'   entry in \code{feature_names}. Default = NULL.
#' 
#' @param feature_names Name of column in \code{rowData} containing names of
#'   continuous features to plot (e.g. gene names). For example, set to
#'   \code{feature_names = "gene_name"} if gene names are stored in a column
#'   named \code{"gene_name"}. This argument is used if \code{annotate} is a
#'   continuous variable. Default = NULL, in which case the row names of the
#'   input object will be used.
#' 
#' @param assay_name Name of \code{assay} in input object containing values to
#'   plot for a continuous variable. Default = "counts".
#' 
#' @param pal Color palette for annotations. Options for discrete values are
#'   "libd_layer_colors", "Okabe-Ito", or any vector of color names or hex
#'   values. For continuous values, provide a vector of length 2 for the low and
#'   high range, e.g. c("gray90", "navy").
#' 
#' @param point_size Point size. Default = 0.3.
#' 
#' @param legend_position Legend position for discrete annotations. Options are
#'   "left", "right", "top", "bottom", and "none". Default = "right".
#' 
#' @param legend_point_size Legend point size for discrete annotations. Default
#'   = 3.
#' 
#' @param show_axes Whether to show axis titles, text, and ticks. Default =
#'   FALSE.
#' 
#' @param y_reverse Whether to reverse y coordinates. This is usually required
#'   for 10x Genomics Visium datasets when using the default coordinate values.
#'   Default = TRUE. Set to FALSE if not needed, e.g. for other platforms.
#' 
#' @param text_by Column name of annotation labels to display over each cluster
#'   of points. This will usually be the same as \code{annotate}. Alternatively,
#'   another column may be used (e.g. with more readable classes or shorter
#'   strings). Only used for discrete \code{annotate}. Default = NULL.
#' 
#' @param text_by_size Text size for annotation labels over each cluster.
#'   Default = 5.
#' 
#' @param text_by_color Color name or hex code for annotation labels. Default =
#'   "black".
#' 
#' 
#' @return Returns a ggplot object, which may be further modified using ggplot
#'   functions.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment reducedDimNames reducedDim
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
#' @importFrom stats median
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes_string geom_point facet_wrap coord_fixed
#'   theme_bw theme element_blank scale_color_viridis_c scale_color_gradientn
#'   scale_color_gradient scale_color_manual ggtitle labs guides scale_y_reverse
#'   aes .data
#' 
#' @export
#' 
#' @author Lukas M. Weber and Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' 
#' # discrete annotations
#' spe <- Visium_humanDLPFC()
#' plotSpots(spe, annotate = "ground_truth")
#' 
#' # continuous annotations
#' spe <- Visium_mouseCoronal()
#' plotSpots(spe, annotate = "Gapdh", feature_names = "gene_name")
#' 
plotSpots <- function(spe, x_coord = NULL, y_coord = NULL, 
                      sample_id = NULL, in_tissue = "in_tissue", 
                      annotate = NULL, feature_names = NULL, 
                      assay_name = "counts", 
                      pal = NULL, point_size = 0.3, 
                      legend_position = "right", 
                      legend_point_size = 3, 
                      show_axes = FALSE, y_reverse = TRUE, 
                      text_by = NULL, text_by_size = 5, 
                      text_by_color = "black") {
  
  # check validity of arguments
  if (!is.null(in_tissue)) {
    stopifnot(is.character(in_tissue))
  }
  
  # get names of continuous features
  if (is.null(feature_names)) {
    feature_nms <- rownames(spe)
  } else {
    stopifnot(ncol(rowData(spe)) > 0)
    stopifnot(feature_names %in% colnames(rowData(spe)))
    feature_nms <- rowData(spe)[, feature_names]
  }
  
  if (!is.null(annotate)) {
    stopifnot(is.character(annotate))
    if (!(annotate %in% c(colnames(colData(spe)), feature_nms))) {
      stop("'annotate' should be either (i) the name of a column in colData ", 
           "or (ii) an entry in either a column named 'feature_names' in ", 
           "rowData or the rownames of the input object")
    }
  }
  
  stopifnot(legend_position %in% c("left", "right", "top", "bottom", "none"))
  
  if (!is.null(sample_id)) {
    stopifnot(sample_id %in% colnames(colData(spe)))
    n_samples <- length(table(colData(spe)[, sample_id]))
  } else {
    n_samples <- NULL
  }
  
  # set up data frame for plotting
  if (is(spe, "SpatialExperiment")) {
    # select default columns of x and y coordinates
    if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
    if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
    df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  } else if (is(spe, "SingleCellExperiment")) {
    if (is.null(x_coord) || is.null(y_coord)) {
      stop("Please provide 'x_coord' and 'y_coord' arguments to specify ", 
           "columns in colData containing x and y coordinates.")
    }
    df <- as.data.frame(colData(spe))
  }
  
  if (!is.null(annotate)) {
    # continuous annotation values
    if (annotate %in% feature_nms) {
      stopifnot(is.character(assay_name))
      ix <- which(feature_nms == annotate)
      df[[annotate]] <- assay(spe, assay_name)[ix, ]
    }
    # discrete annotation values
    if ((annotate %in% colnames(colData(spe))) && 
        (is.character(colData(spe)[, annotate]))) {
      df[[annotate]] <- as.factor(df[[annotate]])
    }
  }
  
  if (!is.null(in_tissue)) {
    df <- df[df[, in_tissue] == 1, ]
  }
  
  # color palettes
  if ((is.null(annotate) || is.numeric(df[[annotate]])) && is.null(pal)) {
    # for continuous values, change NULL to arbitrary color so length(pal) == 1
    pal <- "blue"
  }
  # accepts "libd_layer_colors" and "Okabe-Ito", or arbitrary color palette for NULL
  pal <- .get_pal(pal)
  
  
  # main plot
  
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = annotate)) + 
    geom_point(size = point_size) + 
    coord_fixed() + 
    theme_bw() + 
    theme(legend.position = legend_position, 
          panel.grid = element_blank())
  
  if (!is.null(n_samples) && n_samples > 1) {
    p <- p + facet_wrap(~ sample_id)
  }
  
  
  # additional plot formatting
  
  if (show_axes == FALSE) {
    p <- p + theme(axis.title = element_blank(), 
                   axis.text = element_blank(), 
                   axis.ticks = element_blank())
  }
  
  # color scale
  if (!is.null(annotate)) {
    scaling <- if (is.numeric(df[[annotate]])) {
      # continuous values
      if (length(pal) == 1 && 
          pal %in% c("viridis", "magma", "inferno", "plasma", "cividis", 
                     "rocket", "mako", "turbo")) {
        scale_color_viridis_c(option = pal)
      } else if (length(pal) == 1 && pal == "seuratlike") {
        colors <- colorRampPalette(
          colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
        scale_color_gradientn(colors = colors, limits = range(df[[annotate]]))
      } else {
        scale_color_gradient(low = pal[1], high = pal[2])
      }
    } else if (is.factor(df[[annotate]]) || is.character(df[[annotate]])) {
      # discrete values
      if (is.null(pal)) {
        scale_color_manual(name = annotate, 
                           values = hue_pal()(length(unique(df[[annotate]]))))
      } else if (!is.null(pal)) {
        scale_color_manual(values = pal)
      }
    }
  }
  
  # plot title
  if (!is.null(annotate)) {
    if (is.numeric(df[[annotate]])) {
      # continuous values: display plot title but no legend title
      p <- p + 
        scaling + 
        ggtitle(annotate) + 
        labs(color = NULL) + 
        theme(plot.title = element_text(hjust = 0.5))
    } else if (is.factor(df[[annotate]]) || is.character(df[[annotate]])) {
      # discrete values: display legend title but no plot title
      p <- p + 
        scaling + 
        guides(color = guide_legend(override.aes = list(size = legend_point_size)))
    }
  }
  
  # text annotations
  if (!is.null(annotate)) {
    if (is.factor(df[[annotate]]) || is.character(df[[annotate]])) {
      # add text with the median locations of the 'text_by' vector
      if (!is.null(text_by)) {
        by_text_x <- vapply(
          split(df[[x_coord]], df[[text_by]]), 
          median, 
          FUN.VALUE = 0
        )
        by_text_y <- vapply(
          split(df[[y_coord]], df[[text_by]]), 
          median, 
          FUN.VALUE = 0
        )
        p <- p + 
          geom_text_repel(
            data = data.frame(
              x = by_text_x, 
              y = by_text_y, 
              label = names(by_text_x)
            ), 
            mapping = aes(x = .data$x, y = .data$y, label = .data$label), 
            size = text_by_size, 
            color = text_by_color
          )
      }
    }
  }
  
  # reverse y axis
  if (y_reverse) {
    p <- p + scale_y_reverse()
  }
  
  # return plot
  p
}
