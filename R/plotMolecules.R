#' plotMolecules
#' 
#' Plotting functions for spatially resolved transcriptomics data.
#' 
#' Function to plot molecule-based spatially resolved transcriptomics data
#' stored in a \code{SpatialExperiment} object.
#' 
#' This function generates a plot in spatial coordinates (e.g. x-y coordinates
#' on a tissue slide), for a selected molecule.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param molecule (character) Name of mRNA molecule to plot (assumed to match
#'   one of the row names of \code{rowData}).
#' 
#' @param x_coord (character) Name of column in \code{spatialCoords} containing
#'   x-coordinates of the cell centroids. Default = NULL, which selects the
#'   first column.
#' 
#' @param y_coord (character) Name of column in \code{spatialCoords} containing
#'   y-coordinates of the cell centroids. Default = NULL, which selects the
#'   second column.
#' 
#' @param sample_id (character) Name of column in \code{colData} containing
#'   sample IDs. For datasets with multiple samples, this is used to plot
#'   multiple panels (one per sample) using facetting.
#' 
#' @param palette (character) Color palette, provided as a vector of length 2
#'   for the low and high range. Default = \code{c("gray90", "navy")}.
#' 
#' @param size (numeric) Point size for \code{geom_point()}. Default = 0.3.
#' 
#' 
#' @return Returns a ggplot object. Additional plot elements can be added as
#'   ggplot elements (e.g. title, labels, formatting, etc).
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point
#'   scale_color_gradient coord_fixed ggtitle theme_void
#' 
#' @export
#' 
#' @examples
#' library(STexampleData)
#' spe <- seqFISH_mouseEmbryo()
#' plotMolecules(spe, molecule = "Sox2")
#' 
plotMolecules <- function(spe, 
                          molecule =  NULL, 
                          x_coord = NULL, y_coord = NULL, 
                          sample_id = "sample_id", 
                          assay_name = "counts",
                          palette = c("gray90", "navy"), 
                          size = 0.3) {
  
  if (!is.null(molecule)) stopifnot(is.character(molecule))
  
  if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
  
  n_samples <- length(table(colData(spe)[, sample_id]))
  
  mRNA_counts <- as.numeric(assay(spe, assay_name)[molecule, ])
  stopifnot(length(mRNA_counts) == ncol(spe))
  
  # providing a single value e.g. "navy" will create a vector c("gray95", "navy")
  palette <- .get_pal(palette)
  
  df <- cbind.data.frame(spatialCoords(spe), sum = mRNA_counts)
  
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = "sum")) + 
    geom_point(size = size) +
    coord_fixed() +
    ggtitle(molecule) +
    theme_void()
  
  if (n_samples > 1) {
    p <- p + facet_wrap(~ sample_id)
  }
  
  if(length(palette) == 1 && palette == "viridis"){
    p <- p + scale_colour_viridis_c()
  }else if(length(palette) == 2){
    p <- p + scale_color_gradient(low = palette[1], high = palette[2], trans = "sqrt") 
  }
  
  p
}

