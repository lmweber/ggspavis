#' plotSpots
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot spot-based spatially resolved transcriptomics data stored in
#' a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in spatial coordinates (e.g. x-y coordinates
#' on a tissue slide), along with annotation such as cluster labels or total UMI
#' counts.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} or \code{SingleCellExperiment} object.
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates. Default = NULL, which selects the first column.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates. Default = NULL, which selects the second column.
#' 
#' @param sample_id (character) Name of column in \code{colData} containing
#'   sample IDs. For datasets with multiple samples, this is used to plot
#'   multiple panels (one per sample) using facetting. Default = \code{NULL}.
#' 
#' @param in_tissue (character) Name of column in \code{colData} identifying
#'   spots over tissue, e.g. "in_tissue" for 10x Genomics Visium data. If this
#'   argument is provided, only spots over tissue will be shown. Alternatively,
#'   set to NULL to display all spots. Default = "in_tissue".
#' 
#' @param annotate (character) Name of column in \code{colData} containing
#'   values to annotate spots with colors, e.g. cluster labels (discrete values)
#'   or total UMI counts (continuous values).
#' 
#' @param palette (character) Color palette for annotation. Options for discrete
#'   labels (e.g. cluster labels) are "libd_layer_colors", "Okabe-Ito", or a
#'   vector of color names or hex values. For continuous values (e.g. total UMI
#'   counts), provide a vector of length 2 for the low and high range, e.g.
#'   \code{c("gray90", "navy")}. Default = \code{"libd_layer_colors"}.
#' 
#' @param y_reverse (logical) Whether to reverse y coordinates, which is often
#'   required for 10x Genomics Visium data. Default = TRUE.
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
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point coord_fixed
#'   ggtitle theme_bw theme element_blank scale_y_reverse scale_color_manual
#'   scale_color_gradient
#' 
#' @export
#' 
#' @author Lukas M. Weber with modifications by Yixing E. Dong
#' 
#' @examples
#' library(STexampleData)
#' spe <- Visium_humanDLPFC()
#' plotSpots(spe, annotate = "ground_truth")
#' 
plotSpots <- function(spe, 
                      x_coord = NULL, y_coord = NULL, 
                      sample_id = NULL, 
                      in_tissue = "in_tissue", 
                      trans = "identity",
                      assay = "counts", legend.position = "right", 
                      annotate = NULL, palette = NULL, 
                      y_reverse = TRUE, pt.size = 0.3, show_axis = FALSE,
                      text_by = NULL, text_by_size = 5, text_by_color = "black") {
  
  if (!is.null(in_tissue)) stopifnot(is.character(in_tissue))
  stopifnot(legend.position %in% c("left", "right", "top", "bottom", "none"))
  
  stopifnot(is.character(annotate))
  
  if (!is.null(sample_id)){
    stopifnot(sample_id %in% colnames(colData(spe)))
    n_samples <- length(table(colData(spe)[, sample_id]))
  }else{
    n_samples <- NULL
  }
  
  if(class(spe) == "SingleCellExperiment"){
    if (is.null(x_coord)) {stop("Please specify x_coord name in the SCE colData")}
    if (is.null(y_coord)) {stop("Please specify y_coord name in the SCE colData")}
    
    plt_df <- cbind.data.frame(colData(spe))
  }else if(class(spe) == "SpatialExperiment"){
    if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
    if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
    
    plt_df <- cbind.data.frame(colData(spe), spatialCoords(spe))
  }
  
  if(is.character(plt_df[[annotate]])) plt_df[[annotate]] <- as.factor(plt_df[[annotate]])
  
  if(!annotate %in% c(names(plt_df), rownames(spe))){
    stop("'annotate' should be in rownames(spe) or names(colData(spe))")
  }
  # (optionally) add feature assay data to 'plt_df'
  if(annotate %in% rownames(spe)){
    stopifnot(is.character(assay))
    plt_df[[annotate]] <- assay(spe, assay)[annotate, ]
  }
  
  if (!is.null(in_tissue)) {
    plt_df <- plt_df[plt_df[, in_tissue] == 1, ]
  }
  
  if(is.numeric(plt_df[[annotate]]) & is.null(palette)){
    palette <- "seuratlike" # for continuous feature, turn length(palette) = 0 to length(palette) = 1
  }
  # accepts "libd_layer_colors" and "Okabe-Ito"
  palette <- .get_pal(palette)
  
  p <- ggplot(plt_df, aes_string(x = x_coord, y = y_coord, color = annotate)) + 
    geom_point(size = pt.size) + 
    coord_fixed() + 
    theme_bw() 
  
  if(show_axis == TRUE){
    p <- p + theme(plot.title = element_text(hjust = 0.5),
                   legend.position = legend.position)
  }else{
    p <- p + theme(panel.border = element_blank(),
                   panel.grid = element_blank(),
                   axis.title = element_blank(),
                   axis.text = element_blank(),
                   axis.ticks = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   legend.position = legend.position) 
  }
  
  if (!is.null(n_samples)) {
    if(n_samples > 1){
      p <- p + facet_wrap(~ sample_id)
    }
  }
  
  scale <- if(is.numeric(plt_df[[annotate]])){
    if(length(palette) == 1 && palette == "viridis"){
      scale_color_viridis_c(trans = trans) 
    }else if(length(palette) == 1 && palette == "seuratlike"){
      scale_color_gradientn(colors = colorRampPalette(colors = rev(x = RColorBrewer::brewer.pal(n = 11, name = "Spectral")))(100),
                            trans = trans, limits = c(min(plt_df[[annotate]]), max(plt_df[[annotate]])))
    }else{
      scale_color_gradient(low = palette[1], high = palette[2], trans = trans)
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
        split(plt_df[[x_coord]], plt_df[[text_by]]),
        median,
        FUN.VALUE = 0
      )
      
      by_text_y <- vapply(
        split(plt_df[[y_coord]], plt_df[[text_by]]),
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
  
  if (y_reverse) {
    p <- p + scale_y_reverse()
  }
  
  p
}

