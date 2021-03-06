% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMolecules.R
\name{plotMolecules}
\alias{plotMolecules}
\title{plotMolecules}
\usage{
plotMolecules(
  spe,
  molecule = NULL,
  x_coord = "x",
  y_coord = "y",
  palette = NULL
)
}
\arguments{
\item{spe}{(SpatialExperiment) Input data object.}

\item{molecule}{(character) Name of mRNA molecule to plot (matching to one of
the row names of 'rowData').}

\item{x_coord}{(character) Name of column in 'spatialData' containing
x-coordinates of the cell centroids. Default = 'x'.}

\item{y_coord}{(character) Name of column in 'spatialData' containing
y-coordinates of the cell centroids. Default = 'y'.}

\item{palette}{(character) Color palette for points. Options are a single
color name (e.g. 'red', 'navy', etc), or a vector of length two containing
color names for each end of the scale. Default = 'navy'.}
}
\value{
Returns a ggplot object. Additional plot elements can be added as
  ggplot elements (e.g. title, customized formatting, etc).
}
\description{
Plots for spatially resolved transcriptomics datasets
}
\details{
Function to plot molecule-based spatially resolved transcriptomics data in
spatial (x-y) coordinates.

This function generates a plot showing counts for a given molecule in the x-y
coordinates of the tissue slide.
}
\examples{
# library(ggspavis)
# library(STdata)
# spe <- load_data("seqFISH_mouseEmbryo")
# plotMolecules(spe, molecule = "Sox2")

}
