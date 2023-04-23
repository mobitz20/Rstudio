#https://satijalab.org/seurat/articles/spatial_vignette.html

#remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE, lib = 'D:/Program Files/R-4.2.3/library')
#install.packages("devtools", lib = 'D:/Program Files/R-4.2.3/library')
#devtools::install_github('satijalab/seurat-data')
#install.packages("Rfast2")
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)

InstallData("stxBrain")
?stxBrain
brain <- LoadData("stxBrain", type = "anterior1")

#Data pre-processing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)

#Gene expression visualization
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))

library(ggplot2)
plot <- SpatialFeaturePlot(brain, features = c("Ttr")) + theme(legend.text = element_text(size = 0),
                                                               legend.title = element_text(size = 20), legend.key.size = unit(1, "cm"))
jpeg(filename = "/cloud/project/spatial_vignette_ttr.jpg", height = 700, width = 1200, quality = 50)
print(plot)
dev.off()

#Dimensionality reduction, clustering, and visualizatio
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2

SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,
                                                                                     5, 8)), facet.highlight = TRUE, ncol = 3)
SpatialDimPlot(brain, interactive = TRUE)
SpatialFeaturePlot(brain, features = "Ttr", interactive = TRUE)
LinkedDimPlot(brain)

#Identification of Spatially Variable Features
de_markers <- FindMarkers(brain, ident.1 = 11, ident.2 = 8)
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000],
                                       selection.method = "moransi")

top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 6)
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))


#Subset out anatomical regions
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
# now remove additional cells, use SpatialDimPlots to visualize what to remove
# SpatialDimPlot(cortex,cells.highlight = WhichCells(cortex, expression = image_imagerow > 400
# | image_imagecol < 150))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)

p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
