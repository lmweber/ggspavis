
CD <- cbind.data.frame(spatialCoords(vis), colData(vis))

p <- ggplot(CD, aes(x = array_col, y = array_row, color = as.factor(cluster))) +
  geom_point(size = 0.3) +
  labs(color = NULL) +
  theme_bw() +
  theme(legend.position="right") +
  scale_color_manual(name = "cluster",
                     values = scales::hue_pal()(length(unique(vis[["cluster"]])))) +
  ggtitle("cluster") +
  guides(colour = guide_legend(override.aes = list(size=3)))

p

plotSpots()



library(STexampleData)

spe <- Visium_mouseCoronal()

# color by x coordinate, highlight in-tissue spots
plotVisium(spe, fill = "pxl_col_in_fullres", highlight = "in_tissue")

# subset in-tissue spots
sub <- spe[, as.logical(colData(spe)$in_tissue)]

# color by feature counts, don't include image
rownames(sub) <- make.names(rowData(sub)$gene_name)
plotVisium(sub, annotate = "Gad2", assay = "counts")

plotVisium(sub, annotate = "Xkr4", assay = "counts")



plotMolecules(sub, molecule = "Gad2")


# Test 1, unknown palette
plotVisium(vis, annotate = "cluster_factor")

# Test 2, known palette
library(STexampleData)

spe <- Visium_humanDLPFC()
plotVisium(spe, annotate = "ground_truth", palette = "libd_layer_colors")
plotVisium(spe, annotate = "ground_truth", palette = "libd_layer_colors", highlight = "in_tissue")








