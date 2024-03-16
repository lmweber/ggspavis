#' plotDimRed
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot spot-based spatially resolved transcriptomics data stored in
#' a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in reduced dimension coordinates (PCA or
#' UMAP), along with annotation such as cluster labels or total UMI counts.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} or \code{SingleCellExperiment} object.
#' 
#' @param plot_type (character) Type of reduced dimension plot. Options can be
#'   "UMAP", "PCA", or any self-defined reduced dimension names. Default =
#'   "UMAP".
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values).
#' 
#' @param features (character) Name of column in \code{rowData} containing
#'   continuous feature to plot. Default = "gene_name".
#' 
#' @param assay_name (character) Name of \code{assay} containing continuous
#'   feature to plot.
#' 
#' @param update_colnames (logical) Whether to update column names of
#'   \code{reducedDim} to default values. Default = TRUE.
#' 
#' @param pal (character) Color palette for annotation. Options for discrete
#'   labels are "libd_layer_colors", "Okabe-Ito", or a vector of color names or
#'   hex values. For continuous values, provide a vector of length 2 for the low
#'   and high range, e.g. \code{c("gray90", "navy")}. Default =
#'   \code{"libd_layer_colors"}.
#' 
#' @param point_size (numeric) Point size for \code{geom_point()}. Default = 0.3.
#' 
#' @param text_by (character) Column name of the annotation to apply on top of 
#' each cluster. Usually should put it the same as `annotate = `. unless you have 
#' another intended `text_by` column, e.g. with more readable classes or shorter 
#' strings. Only used for categorical `annotate = `. Default = \code{NULL}.
#' 
#' @param text_by_size (numerical) Text size hovering over each cluster. Default =
#'  \code{5}.
#' 
#' @param text_by_color (character) Color string or hex code. Default = \code{"black"}.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SingleCellExperiment reducedDimNames reducedDim colData rowData
#' @importFrom SummarizedExperiment assay
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom scales hue_pal
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab theme_bw theme
#'   element_blank scale_color_viridis_c scale_color_gradientn
#'   scale_color_gradient scale_color_manual ggtitle labs guides
#' 
#' @export
#' 
#' @author Lukas M. Weber with modifications by Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' 
#' # use small subset of data for this example
#' # for longer examples see our online book OSTA
#' spe <- spe[, colData(spe)$in_tissue == 1]
#' set.seed(100)
#' n <- 200
#' spe <- spe[, sample(seq_len(ncol(spe)), n)]
#' 
#' # calculate log-transformed normalized counts
#' library(scran)
#' set.seed(100)
#' qclus <- quickCluster(spe)
#' spe <- computeSumFactors(spe, cluster = qclus)
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
#' set.seed(100)
#' spe <- runPCA(spe, subset_row = top_hvgs)
#' set.seed(100)
#' spe <- runUMAP(spe, dimred = "PCA")
#' colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)
#' 
#' # generate plot
#' plotDimRed(spe, type = "UMAP", annotate = "ground_truth")
#' 
plotDimRed <- function(spe, plot_type = c("UMAP", "PCA"), 
                       annotate = NULL, features = "gene_name", 
                       assay_name = "counts", 
                       update_colnames = TRUE, 
                       palette = NULL, point_size = 0.3, 
                       text_by = NULL, text_by_size = 5, 
                       text_by_color = "black") {
  
  # check validity of arguments
  plot_type <- match.arg(plot_type)
  stopifnot(plot_type %in% reducedDimNames(spe))
  
  stopifnot(is.character(annotate))
  if (!(annotate %in% c(colnames(colData(spe)), rowData(spe)[, features]))) {
    stop("'annotate' should be the name of a column in colData or an entry in ", 
         "the column 'features' in rowData")
  }
  
  # update colnames of reducedDims for plotting
  if (update_colnames) {
    colnames(reducedDim(spe, plot_type)) <- 
      paste0(plot_type, "_", seq_len(ncol(reducedDim(spe, plot_type))))
  }
  
  # data for plotting
  df <- cbind.data.frame(colData(spe), reducedDim(spe, plot_type))
  
  # continuous values
  if (annotate %in% rowData(spe)[, features]) {
    stopifnot(is.character(assay_name))
    ix <- which(rowData(spe)[, features] == annotate)
    df[[annotate]] <- assay(spe, assay_name)[ix, ]
  }
  if ((annotate %in% colnames(colData(spe))) && 
      (is.character(colData(spe)[, annotate]))) {
    df[[annotate]] <- as.factor(df[[annotate]])
  }
  
  # color palettes
  if (is.numeric(df[[annotate]]) && is.null(palette)) {
    # for continuous feature, change NULL to an arbitrary color so length(pal) == 1
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
  scaling <- if (is.numeric(df[[annotate]])) {
    # continuous values
    if (length(pal) == 1 && 
        palette %in% c("viridis", "magma", "inferno", "plasma", "cividis", 
                       "rocket", "mako", "turbo")) {
      scale_color_viridis_c(option = pal)
    } else if (length(pal) == 1 && palette == "seuratlike") {
      colors <- colorRampPalette(
        colors = rev(x = brewer.pal(n = 11, name = "Spectral")))(100)
      scale_color_gradientn(
        colors = colorRampPalette(colors = colors), 
        limits = range(df[[annotate]]))
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
  
  # plot title
  if (is.numeric(df[[annotate]])) {
    # continuous values: display plot title but no legend title
    p <- p + 
      scaling + 
      ggtitle(annotate) + 
      labs(color = NULL)
  } else if (is.factor(df[[annotate]]) | is.character(df[[annotate]])) {
    # discrete values: display legend title but no plot title
    p <- p + 
      scaling + 
      guides(color = guide_legend(override.aes = list(size = 3)))
  }
  
  # text annotations
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
  
  # return plot
  p
}
