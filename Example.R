library(SpatialExperiment)
library(ggspavis)
library(grid)
library(scater)

path <- "/home/estelladong/Desktop/SampleData/archive/Xenium_Preprint_Data/Visium/outs"
vis <- SpatialExperiment::read10xVisium(path)

# add some values in 'colData' to annotate spots
colData(vis)$sum <- colSums(counts(vis))
vis <- scater::logNormCounts(vis)

plotVisium(vis, fill = "sum", palette = "seuratlike")
ggspavis::plotVisium(vis, fill = "sum")

plotMolecules(vis, molecule = "ENSG00000009307", assay = "logcounts")

plotMolecules(vis, molecule = "ENSG00000009307", assay = "counts", palette = "viridis")
plotMolecules(vis, molecule = "ENSG00000009307", assay = "counts", palette = "seuratlike")
plotMolecules(vis, molecule = "ENSG00000009307", assay = "logcounts", palette = "seuratlike")


plotMolecules(vis, molecule = "ENSG00000009307", assay = "counts")
ggspavis::plotMolecules(vis, molecule = "ENSG00000009307")

library(dplyr)
as.data.frame(rowData(vis)) %>%
  filter(symbol == "CSDE1")


vis$cluster <- sample(1:10, ncol(vis), replace=T)
plotVisium(vis, fill = "cluster", palette = "seuratlike")

# plotVisium() is for total count continuous 
# plotVisium() is also for cluster categorical (any levels)
# plotMolecules() is for genes expression continuous -> redundant!!


# plot

vis$cluster_factor <- factor(vis$cluster)
plotVisium(vis, fill = "cluster_factor")


plotSpots(vis, annotate = "cluster")
plotSpots(vis, annotate = "cluster_factor")

plotVisium(vis, fill = "cluster_factor")



1. Remove plotMolecules()
plot Visium is for with HE on Visium only, plot Spot is for without HE and any technology, xe, cos, mer, vis, etc.
2. plotVisium add option for any number of factor cluster with no default palette
2. plotVisium add option for two more palettes ("viridis" and "seuratlike") on cont var (sum or any gene) with default "seuratlike"
# Maybe need to rename seuratlike as rainbow, and put as default 
3. same update for plotSpots and plotDimRed

3. unify the feature name (plotVisium::fill, plotSpots::annotate, plotDimRed::annotate) as just annotate
4. Unify df to all plt_df inside the functions






