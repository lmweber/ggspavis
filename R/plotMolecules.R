#' plotMolecules
#' 
#' Plotting functions for spatial transcriptomics data.
#' 
#' Function to create spot plot for molecule-based datasets, showing spatial
#' locations in x-y coordinates with optional annotations such as expression of
#' a gene.
#' 
#' 
#' @param spe (SpatialExperiment) Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param molecule Name of mRNA molecule to plot (assumed to match one of the
#'   row names of \code{rowData}).
#' 
#' @param x_coord Name of column in \code{spatialCoords} containing x
#'   coordinates. Default = NULL, which selects the first column of
#'   \code{spatialCoords}.
#' 
#' @param y_coord Name of column in \code{spatialCoords} containing y
#'   coordinates. Default = NULL, which selects the second column of
#'   \code{spatialCoords}.
#' 
#' @param sample_id Name of column in \code{colData} containing sample IDs. This
#'   argument is only required for datasets containing multiple samples (tissue
#'   sections). If provided, samples will be shown in multiple panels using
#'   facetting. Default = NULL.
#' 
#' @param pal Color palette, provided as a vector of length 2 for the low and
#'   high range. Default = c("gray90", "navy").
#' 
#' @param point_size Point size. Default = 0.3.
#' 
#' 
#' @return Returns a ggplot object, which may be further modified using ggplot
#'   functions.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment counts
#' @importFrom ggplot2 ggplot aes_string facet_wrap geom_point
#'   scale_color_gradient coord_fixed ggtitle theme_void
#' 
#' @export
#' 
#' @author Lukas M. Weber
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
                          pal = c("gray90", "navy"), 
                          point_size = 0.3) {
  
  if (!is.null(molecule)) stopifnot(is.character(molecule))
  
  if (is.null(x_coord)) x_coord <- colnames(spatialCoords(spe))[1]
  if (is.null(y_coord)) y_coord <- colnames(spatialCoords(spe))[2]
  
  n_samples <- length(table(colData(spe)[, sample_id]))
  
  mRNA_counts <- as.numeric(counts(spe)[molecule, ])
  stopifnot(length(mRNA_counts) == ncol(spe))
  
  # providing a single value e.g. "navy" will create a vector c("gray95", "navy")
  pal <- .get_pal(pal)
  
  df <- cbind.data.frame(spatialCoords(spe), sum = mRNA_counts)
  
  p <- ggplot(df, aes_string(x = x_coord, y = y_coord, color = "sum")) + 
    geom_point(size = size) + 
    scale_color_gradient(low = palette[1], high = palette[2], trans = "sqrt") + 
    coord_fixed() + 
    ggtitle(molecule) + 
    theme_void()
  
  if (n_samples > 1) {
    p <- p + facet_wrap(~ sample_id)
  }
  
  p
}
