#' plotVisium
#' 
#' Plots for spatially resolved transcriptomics data from the 10x Genomics
#' Visium platform
#' 
#' Function to plot spatially resolved transcriptomics data from the spot-based
#' 10x Genomics Visium platform.
#' 
#' This function generates a plot for spot-based spatially resolved
#' transcriptomics data from the 10x Genomics Visium platform, with several
#' options available to adjust the plot type and style.
#' 
#' 
#' @param spe (SpatialExperiment) Input data object.
#' 
#' @param spots (logical) Whether to display spots (spatial barcodes) as points.
#'   Default = TRUE.
#' 
#' @param fill (character) Column in 'colData' to use to fill points by color.
#'   If 'fill' contains a numeric column (e.g. total UMI counts), a continuous
#'   color scale will be used. If 'fill' contains a factor (e.g. cluster
#'   labels), a discrete color scale will be used. Default = NULL.
#' 
#' @param highlight (character) Column in 'spatialData' to use to highlight
#'   points by outlining them. For example, 'in_tissue' will highlight spots
#'   overlapping with tissue. Default = 'in_tissue'. Set to NULL to disable.
#' 
#' @param facets (character) Column in 'colData' to use to facet plots, i.e.
#'   show multiple panels of plots. Default = 'sample_id'. Set to NULL to
#'   disable.
#' 
#' @param image (logical) Whether to show histology image as background. Default
#'   = TRUE.
#' 
#' @param x_coord (character) Column in 'spatialData' containing x-coordinates.
#'   Default = 'x'.
#' 
#' @param y_coord (character) Column in 'spatialData' containing y-coordinates.
#'   Default = 'y'.
#' 
#' @param flip_xy_Visium (logical) Whether to flip x and y coordinates and
#'   reverse y scale to match orientation of histology images. This is sometimes
#'   required for Visium data, depending on the orientation of the raw data.
#' 
#' @param palette (character) Color palette for points. Options for discrete
#'   labels are 'libd_layer_colors', 'Okabe-Ito', or a custom vector of hex
#'   color codes. Options for continuous values are 'viridis', a single color
#'   name (e.g. 'red', 'navy', etc), or a vector of length two containing color
#'   names for each end of the scale. Default = 'libd_layer_colors' for discrete
#'   data, and 'viridis' for continuous data.
#' 
#' @param sample_ids (character) Samples to show, if multiple samples are
#'   available. Default = NULL (show all samples).
#' 
#' @param image_ids (character) Images to show, if multiple images are
#'   available. Default = NULL (show all images).
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, customized formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialData imgData 'imgData<-' scaleFactors
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual scale_fill_gradient
#'   scale_fill_viridis_c scale_color_identity facet_wrap guides guide_colorbar
#'   guide_legend theme_void element_text margin unit layer
#' @importFrom methods is as
#' @importFrom stats setNames
#' 
#' @export
#' 
#' @author Helena L. Crowell and Lukas M. Weber
#' 
#' @examples
#' # library(ggspavis)
#' # library(STdata)
#' # spe <- load_data("Visium_mouseCoronal")
#' # plotVisium(spe)
#' 
plotVisium <- function(spe, 
                       spots = TRUE, fill = NULL, highlight = "in_tissue", 
                       facets = "sample_id", image = TRUE, 
                       x_coord = "x", y_coord = "y", 
                       flip_xy_Visium = FALSE, 
                       sample_ids = NULL, image_ids = NULL, 
                       palette = NULL) {
  
  # check input type
  stopifnot(is(spe, "SpatialExperiment"))
  
  # set up color palette
  if (!is.null(fill)) {
    if (is.null(palette) && is.factor(colData(spe)[[fill]])) {
      palette <- "libd_layer_colors" }
    if (is.null(palette) && is.numeric(colData(spe)[[fill]])) {
      palette <- "viridis" } }
  # if length(palette) > 1, use palette as provided (either multiple colors for
  # discrete labels, or length 2 for continuous gradient)
  if (!is.null(palette) & length(palette) == 1) {
    if (palette == "libd_layer_colors") {
      palette <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
    } else if (palette == "Okabe-Ito") {
      palette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    } else if (palette == "viridis") {
      # if viridis, use 'scale_fill_viridis_c' instead of hex codes or color names
    } else {
      # if providing a single color name (e.g. 'navy' or 'red'), combine with
      # 'gray95' for continuous color scale
      palette <- c("gray95", palette) } }
  
  # select samples
  if (is.null(sample_ids)) {
    # default to using all samples
    sample_ids <- unique(colData(spe)$sample_id)
  }
  spe <- spe[, colData(spe)$sample_id %in% sample_ids]
  
  # select images
  if (is.null(image_ids)) {
    # default to first available image for each sample
    idx <- SpatialExperiment:::.get_img_idx(spe, TRUE, NULL)
  } else {
    if (length(image_ids) == 1) {
      idx <- SpatialExperiment:::.get_img_idx(spe, TRUE, image_ids)
    } else {
      stopifnot(length(image_ids) == length(sample_ids))
      idx <- mapply(s = sample_ids, i = image_ids, 
                    function(s, i) SpatialExperiment:::.get_img_idx(spe, s, i))
    }
  }
  # subset selected images
  imgData(spe) <- imgData(spe)[idx, ]
  
  # data frame for plotting
  df_plot <- as.data.frame(spatialData(spe))
  df_plot <- cbind(df_plot, sample_id = colData(spe)$sample_id)
  if (!is.null(fill)) {
    df_plot <- cbind(df_plot, colData(spe)[, fill, drop = FALSE])
  } else {
    fill <- "value"
    df_plot[[fill]] <- factor(1)
  }
  
  # scale spatial coordinates
  sfs <- setNames(scaleFactors(spe, sample_ids, image_ids), sample_ids)
  for (i in seq_along(sample_ids)) {
    ix <- df_plot$sample_id == sample_ids[i]
    df_plot[ix, c(x_coord, y_coord)] <- df_plot[ix, c(x_coord, y_coord)] * sfs
    if (flip_xy_Visium) {
      # if required: flip x and y coordinates and reverse y scale to match
      # orientation of images (sometimes required for Visium data)
      x_coord_tmp <- df_plot[ix, y_coord]
      y_coord_tmp <- df_plot[ix, x_coord]
      y_coord_tmp <- (-1 * y_coord_tmp) + min(y_coord_tmp) + max(y_coord_tmp)
      df_plot[ix, x_coord] <- x_coord_tmp
      df_plot[ix, y_coord] <- y_coord_tmp
    }
  }
  
  # construct image layers
  # note: images could also be plotted using 'annotation_custom()', however this
  # does not allow for faceting, so we instead construct a separate image layer
  # for each sample
  if (image) {
    # split imgData by sample
    dfs <- split(imgData(spe), imgData(spe)$sample_id)
    # construct separate image layer for each sample
    images <- lapply(dfs, function(.) layer(
      data = as.data.frame(.), 
      inherit.aes = FALSE, 
      stat = "identity", 
      position = "identity", 
      geom = ggplot2::GeomCustomAnn, 
      params = list(
        #grob = imgGrob(.$data[[1]]), 
        xmin = 0, xmax = .$width, 
        ymin = 0, ymax = .$height)))
  } else {
    images <- NULL
  }
  
  # construct points and highlights
  if (spots) {
    # check whether 'fill' is continuous (numeric) or discrete (factor)
    guide <- ifelse(is.numeric(df_plot[[fill]]), guide_colorbar, guide_legend)
    points <- list(
      guides(fill = guide(
        title = fill, order = 1, override.aes = list(col = NA, size = 3))), 
      geom_point(shape = 21, size = 1, stroke = 0.25, alpha = 0.5))
    if (!is.null(highlight)) {
      df_plot$highlight <- as.factor(df_plot[[highlight]])
      highlights <- list(
        scale_color_manual(highlight, values = c("gray50", "black")), 
        guides(col = guide_legend(override.aes = list(
          size = 2, stroke = 1, col = c("gray", "black")[seq_along(unique(df_plot$highlight))]))))
    } else {
      df_plot$highlight <- "transparent"
      highlights <- scale_color_identity()
    }
  } else {
    # this is required, else the image layer doesn't show
    points <- geom_point(col = "transparent")
    highlights <- NULL
  }
  
  # construct facets
  if (!is.null(facets)) {
    facets <- facet_wrap(facets)
  } else {
    facets <- NULL
  }
  
  # create plot
  p <- 
    ggplot(df_plot, aes_string(x_coord, y_coord, fill = fill, col = "highlight")) + 
      images + points + highlights + facets + 
      coord_fixed() + 
      theme_void() + 
      theme(legend.key.size = unit(0.5, "lines"), 
            strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines")))
  
  # colors
  if (!is.null(fill)) {
    if (is.factor(colData(spe)[[fill]])) {
      p <- p + scale_fill_manual(values = palette)
    }
    if (is.numeric(colData(spe)[[fill]])) {
      if (length(palette) == 1 && palette == "viridis") {
        p <- p + scale_fill_viridis_c(trans = "log10")
      } else {
        p <- p + scale_fill_gradient(low = palette[1], high = palette[2], trans = "log10")
      }
    }
  }
  
  # display plot
  p
}

