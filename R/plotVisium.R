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
#' @param assay (character) Name of assay data to use 
#'   when \code{fill} is in \code{rownames(spe)}.
#'   Should be one of \code{assayNames(spe)}.
#'   
#' @param trans Transformation to apply for continuous scales.
#'   Ignored unless \code{fill} is numeric, e.g. feature expression.
#'   (See \code{\link{ggplot2}{continuous_scale}} for valid options)
#' 
#' @param x_coord (character) Column in 'spatialData' containing x-coordinates.
#'   Default = 'x'.
#' 
#' @param y_coord (character) Column in 'spatialData' containing y-coordinates.
#'   Default = 'y'.
#' 
#' @param flip_xy (logical) Whether to flip x and y coordinates and
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
#' @importFrom SpatialExperiment spatialData spatialCoordsNames imgData
#'   'imgData<-' imgRaster scaleFactors
#' @importFrom SingleCellExperiment colData
#' @importFrom ggplot2 ggplot aes_string scale_fill_manual scale_fill_gradient
#'   scale_fill_viridis_c scale_color_identity facet_wrap guides guide_colorbar
#'   guide_legend theme_void element_text margin unit layer
#' @importFrom grid rasterGrob
#' @importFrom methods is as
#' @importFrom stats setNames
#' 
#' @export
#' 
#' @author Helena L. Crowell and Lukas M. Weber
#' 
#' @examples
#' library(ggspavis)
#' library(STexampleData)
#' 
#' spe <- load_data("Visium_mouseCoronal")
#' 
#' # color by x coordinate, highlight in-tissue spots
#' plotVisium(spe, fill = "x", highlight = "in_tissue")
#' 
#' # subset in-tissue spots
#' sub <- spe[, inTissue(spe)]
#' 
#' # color by feature counts, don't include image
#' rownames(sub) <- make.names(rowData(sub)$gene_name)
#' plotVisium(sub, fill = "Gad2", assay = "counts", image = FALSE)
#' 
plotVisium <- function(spe, 
                       spots = TRUE, fill = NULL, highlight = NULL, 
                       facets = "sample_id", image = TRUE, 
                       assay = "logcounts", trans = "identity",
                       x_coord = "x", y_coord = "y", flip_xy = FALSE, 
                       sample_ids = NULL, image_ids = NULL, palette = NULL) {
  
  # check validity of input arguments
  stopifnot(
    is(spe, "SpatialExperiment"),
    is.logical(spots), length(spots) == 1,
    is.logical(image), length(image) == 1,
    is.logical(flip_xy), length(flip_xy) == 1,
    is.character(x_coord), length(x_coord) == 1,
    is.character(y_coord), length(y_coord) == 1,
    c(x_coord, y_coord) %in% spatialCoordsNames(spe))
  
  # set up data for plotting
  plt_df <- data.frame(spatialData(spe), colData(spe))
  if (!is.null(fill)) {
    # check validity of 'fill' argument
    stopifnot(is.character(fill), length(fill) == 1)
    if (!fill %in% c(names(plt_df), rownames(spe))) {
      stop("'fill' should be in rownames(spe) or",
        " names(colData(spe)/spatialData(spe))")
    }
    # (optionally) add feature assay data to 'plt_df'
    if (fill %in% rownames(spe)) {
      stopifnot(
          is.character(assay), 
          length(grep(assay, assayNames(spe))) == 1)
      plt_df[[fill]] <- assay(spe, assay)[fill, ]
    }
    # get color palette
    palette <- .get_pal(palette, plt_df[[fill]])
  } else {
    fill <- "foo"
    plt_df[[fill]] <- "black"
  }

  if (is.null(sample_ids)) {
    # default to using all samples
    sample_ids <- unique(spe$sample_id)
  } else {
    # subset specified samples
    spe <- spe[, spe$sample_id %in% sample_ids]
  }

  # subset selected images
  img_df <- .sub_imgData(spe, sample_ids, image_ids)
  rownames(img_df) <- img_df$sample_id
  
  # scale spatial coordinates
  for (s in sample_ids) {
    ix <- plt_df$sample_id == s
    xy <- c(x_coord, y_coord)
    sf <- img_df[s, "scaleFactor"]
    plt_df[ix, xy] <- sf * plt_df[ix, xy]
    # if required: flip x and y coordinates and reverse y scale to match
    # orientation of images (sometimes required for Visium data)
    if (flip_xy) plt_df <- .flip_xy(plt_df, x_coord, y_coord)
  }
  
  # construct image layers
  # note: images could also be plotted using 'annotation_custom()', 
  # however, this does not allow for faceting, so we instead
  # construct a separate image layer for each sample
  if (image) {
    images <- lapply(sample_ids, function(s) {
      spi <- img_df[s, "data"]
      img <- imgRaster(spi[[1]])
      layer(
        data = data.frame(sample_id = s),
        inherit.aes = FALSE, 
        stat = "identity", 
        position = "identity", 
        geom = ggplot2::GeomCustomAnn, 
        params = list(
          grob = rasterGrob(img),
          xmin = 0, xmax = ncol(img),
          ymin = 0, ymax = nrow(img))
      )
    })
    xlim <- c(0, ncol(img))
    ylim <- c(0, nrow(img))
  } else images <- xlim <- ylim <- NULL
  
  # construct points and highlights
  if (spots) {
    # check whether 'fill' is continuous (numeric) or discrete (factor)
    guide <- ifelse(is.numeric(plt_df[[fill]]), guide_colorbar, guide_legend)
    points <- list(
      guides(fill = guide(
        title = fill, order = 1, override.aes = list(col = NA, size = 3))), 
      geom_point(shape = 21, size = 1, stroke = 0.25, alpha = 0.5))
    if (!is.null(highlight)) {
      plt_df$highlight <- as.factor(plt_df[[highlight]])
      highlights <- list(
        scale_color_manual(highlight, values = c("gray50", "black")), 
        guides(col = guide_legend(override.aes = list(
          size = 2, stroke = 1, col = c("gray", "black")[
            seq_along(unique(plt_df$highlight))]))))
    } else {
      plt_df$highlight <- "transparent"
      highlights <- scale_color_identity()
    }
  } else {
    # this is required, else the image layer doesn't show
    points <- geom_point(col = "transparent")
    highlights <- NULL
  }
  
  # color scale
  scale <- if (fill != "foo") {
    if (is.numeric(plt_df[[fill]])) {
      if (length(palette) == 1 && palette == "viridis") {
        scale_fill_viridis_c(trans = trans)
      } else {
        scale_fill_gradient(low = palette[1], high = palette[2], trans = trans)
      }
    } else {
      scale_fill_manual(values = palette)
    }
  } else scale_fill_identity()

  # display plot
  ggplot(plt_df, 
    aes_string(x_coord, y_coord, fill = fill, col = "highlight")) + 
    images + points + highlights + scale + 
    coord_fixed(xlim = xlim, ylim = ylim) + 
    theme_void() + theme(
      legend.key.size = unit(0.5, "lines"), 
      strip.text = element_text(margin = margin(0, 0, 0.5, 0, "lines"))) +
    if (!is.null(facets)) facet_wrap(facets)
}
