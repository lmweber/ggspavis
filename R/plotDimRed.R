#' plotDimRed
#' 
#' Plotting functions for spatial transcriptomics data.
#' 
#' Function to create reduced dimension plot (e.g. PCA or UMAP) with additional
#' optional annotations such as cluster labels, expression of a gene, or quality
#' control metrics.
#' 
#' 
#' @param spe Input data, assumed to be a \code{SpatialExperiment} or
#'   \code{SingleCellExperiment} object.
#' 
#' @param plot_type Type of reduced dimension plot. Possible options are "UMAP",
#'   "PCA", or any other set of reduced dimensions stored in the input object.
#'   Default = "UMAP".
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
#' @param update_dimnames Whether to update column names of \code{reducedDims}
#'   to default values for plotting. Default = TRUE.
#' 
#' @param pal Color palette for annotations. Options for discrete values are
#'   "libd_layer_colors", "Okabe-Ito", or any vector of color names or hex
#'   values. For continuous values, provide a vector of length 2 for the low and
#'   high range, e.g. c("gray90", "navy").
#' 
#' @param point_size Point size. Default = 0.3.
#' 
#' @param legend_point_size Legend point size for discrete annotations. Default
#'   = 3.
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
#' @importFrom SingleCellExperiment reducedDimNames reducedDim 'reducedDim<-'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
#' @importFrom stats median
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab theme_bw theme
#'   element_blank scale_color_viridis_c scale_color_gradientn
#'   scale_color_gradient scale_color_manual ggtitle labs guides aes .data
#' 
#' @export
#' 
#' @author Lukas M. Weber and Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' 
#' # select spots over tissue
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' 
#' # use small subset of data for this example
#' n <- 200
#' set.seed(123)
#' spe <- spe[, sample(seq_len(ncol(spe)), n)]
#' 
#' # calculate logcounts
#' library(scran)
#' spe <- logNormCounts(spe)
#' 
#' # identify top highly variable genes (HVGs)
#' is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
#' spe <- spe[!is_mito, ]
#' dec <- modelGeneVar(spe)
#' top_hvgs <- getTopHVGs(dec, prop = 0.1)
#' 
#' # run dimensionality reduction
#' library(scater)
#' set.seed(123)
#' spe <- runPCA(spe, subset_row = top_hvgs)
#' set.seed(123)
#' spe <- runUMAP(spe, dimred = "PCA")
#' colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
#' 
#' # generate plot
#' plotDimRed(spe, plot_type = "UMAP", annotate = "ground_truth")
#' 
plotDimRed <- function(spe, plot_type = c("UMAP", "PCA"), 
                       annotate = NULL, feature_names = NULL, 
                       assay_name = "counts", 
                       update_dimnames = TRUE, 
                       pal = NULL, point_size = 0.3, 
                       legend_point_size = 3, 
                       text_by = NULL, text_by_size = 5, 
                       text_by_color = "black") {
  
  # check validity of arguments
  stopifnot(length(plot_type) == 1)
  stopifnot(plot_type %in% reducedDimNames(spe))
  
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
  
  # update colnames of reducedDims for plotting
  if (update_dimnames) {
    colnames(reducedDim(spe, plot_type)) <- 
      paste0(plot_type, "_", seq_len(ncol(reducedDim(spe, plot_type))))
  }
  
  # data frame for plotting
  df <- cbind.data.frame(colData(spe), reducedDim(spe, plot_type))
  
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
  
  # color palettes
  if ((is.null(annotate) || is.numeric(df[[annotate]])) && is.null(pal)) {
    # for continuous values, change NULL to arbitrary color so length(pal) == 1
    pal <- "blue"
  }
  # accepts "libd_layer_colors" and "Okabe-Ito", or arbitrary color palette for NULL
  pal <- .get_pal(pal)
  
  dim_labels <- colnames(reducedDim(spe, plot_type))[1:2]
  x_label <- dim_labels[1]
  y_label <- dim_labels[2]
  
  
  # main plot
  
  p <- ggplot(df, aes_string(x = x_label, y = y_label, color = annotate)) + 
    geom_point(size = point_size) + 
    xlab(x_label) + 
    ylab(y_label) + 
    theme_bw() + 
    theme(panel.grid = element_blank())
  
  
  # additional plot formatting
  
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
    } else if (is.factor(df[[annotate]]) | is.character(df[[annotate]])) {
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
    } else if (is.factor(df[[annotate]]) | is.character(df[[annotate]])) {
      # discrete values: display legend title but no plot title
      p <- p + 
        scaling + 
        guides(color = guide_legend(override.aes = list(size = legend_point_size)))
    }
  }
  
  # text annotations
  if (!is.null(annotate)) {
    if (is.factor(df[[annotate]]) | is.character(df[[annotate]])) {
      # add text with the median locations of the 'text_by' vector
      if (!is.null(text_by)) {
        by_text_x <- vapply(
          split(df[[x_label]], df[[text_by]]), 
          median, 
          FUN.VALUE = 0
        )
        by_text_y <- vapply(
          split(df[[y_label]], df[[text_by]]), 
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
  
  # return plot
  p
}
