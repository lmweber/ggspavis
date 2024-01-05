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
#' @param type (character) Type of reduced dimension plot. Options can be "UMAP",
#' "PCA", or any self-defined reduced dimension names. Default = "UMAP".
#' 
#' @param x_axis (character) Name of column in \code{reducedDim} containing
#'   x-coordinates. Default = "UMAP1" or "PC1", depending on plot type.
#' 
#' @param y_axis (character) Name of column in \code{reducedDim} containing
#'   y-coordinates. Default = "UMAP2" or "PC2", depending on plot type.
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values).
#' 
#' @param palette (character) Color palette for annotation. Options for discrete
#'   labels are "libd_layer_colors", "Okabe-Ito", or a vector of color names or
#'   hex values. For continuous values, provide a vector of length 2 for the low
#'   and high range, e.g. \code{c("gray90", "navy")}. Default =
#'   \code{"libd_layer_colors"}.
#' 
#' @param pt.size (numeric) Point size for \code{geom_point()}. Default = 0.3.
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
#' @importFrom SingleCellExperiment colData reducedDim
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes_string geom_point xlab ylab ggtitle theme_bw
#'   theme element_blank scale_color_manual scale_color_gradient
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
plotDimRed <- function(spe, 
                       type = c("UMAP", "PCA"), 
                       assay = "counts",
                       annotate = NULL, palette = NULL, 
                       pt.size = 0.3,
                       text_by = NULL, text_by_size = 5, text_by_color = "black") {
  
  # check validity of input arguments
  if(length(type) != 1){
    stop("Please specify reduced dimension type.")
  }
  # check validity of input arguments
  if(!type %in% reducedDimNames(spe)){
    stop("Reduced dimension 'type' does not exist in reducedDimNames(spe).")
  }
  
  stopifnot(is.character(annotate))
  # no matter colnames(reducedDim()) null or not, reorganize for plotting
  colnames(reducedDim(spe, type)) <- paste0(type, "_", 1:ncol(reducedDim(spe, type)))
  
  plt_df <- cbind.data.frame(colData(spe), reducedDim(spe, type))
  
  if(!annotate %in% c(names(plt_df), rownames(spe))){
    stop("'annotate' should be in rownames(spe) or names(colData(spe)).")
  }
  # (optionally) add feature assay data to 'plt_df'
  if(annotate %in% rownames(spe)){
    stopifnot(is.character(assay))
    plt_df[[annotate]] <- assay(spe, assay)[annotate, ]
  }
  
  if(is.numeric(plt_df[[annotate]]) & is.null(palette)){
    palette <- "blue" # for continuous feature, turn length(palette) = 0 to length(palette) = 1
  }
  # accepts "libd_layer_colors" and "Okabe-Ito", or arbitrary color palette for NULL 
  palette <- .get_pal(palette)
  
  p <- ggplot(plt_df, aes_string(x = paste0(type, "_1"), 
                                 y = paste0(type, "_2"), color = annotate)) + 
    geom_point(size = pt.size) + 
    xlab(paste0(type, "_1")) + 
    ylab(paste0(type, "_2")) + 
    # ggtitle("Reduced dimensions") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          #axis.text = element_blank(), 
          #axis.ticks = element_blank(),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5),
          axis.line.x = element_line(size = 0.5, linetype = "solid", colour = "black"),
          axis.line.y = element_line(size = 0.5, linetype = "solid", colour = "black")
          )

  
  scale <- if(is.numeric(plt_df[[annotate]])){
    if(length(palette) == 1 && palette %in% c("viridis", "magma", "inferno", "plasma",
                                              "cividis", "rocket", "mako", "turbo")){
      scale_color_viridis_c(option = palette)
    }else if(length(palette) == 1 && palette == "seuratlike"){
      scale_color_gradientn(colors = colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
                           limits = c(min(plt_df[[annotate]]), max(plt_df[[annotate]])))
    }else{
      scale_color_gradient(low = palette[1], high = palette[2])
    }
  }else if(is.factor(plt_df[[annotate]])){
    if(is.null(palette)){ # for categorical feature, automate palette
      scale_color_manual(name = annotate,
                        values = scales::hue_pal()(length(unique(plt_df[[annotate]])))) 
    }else if(!is.null(palette)){
      scale_color_manual(values = palette)
    }
  }
  
  
  if(is.numeric(plt_df[[annotate]])){ # continuous display plot title but no legend title
    p <- p + scale + ggtitle(annotate) + labs(color = NULL) 
  }else if(is.factor(plt_df[[annotate]])){ # categorical disply legend title but no plot title
    p <- p + scale + 
      guides(colour = guide_legend(override.aes = list(size = 3)))
  }
  
  
  if(is.factor(plt_df[[annotate]])){
    # Adding text with the median locations of the 'text_by' vector.
    if(!is.null(text_by)){
      by_text_x <- vapply(
        split(plt_df[[paste0(type, "_1")]], plt_df[[text_by]]),
        median,
        FUN.VALUE = 0
      )
      
      by_text_y <- vapply(
        split(plt_df[[paste0(type, "_2")]], plt_df[[text_by]]),
        median,
        FUN.VALUE = 0
      )
      
      p <- p +
        ggrepel::geom_text_repel(
          data = data.frame(
            x = by_text_x, y = by_text_y, label = names(by_text_x)
          ),
          mapping = aes(x = .data$x, y = .data$y, label = .data$label),
          size = text_by_size, colour = text_by_color
        )
    }
  }

  
  p
}

