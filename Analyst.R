#project 4 analyst 

#load packages 
library(tidyverse)
library(Seurat)
library(dplyr)
library(patchwork)


# load in the data 
cells <- readRDS("/projectnb/bf528/project_4_scrnaseq/GSM2230760_seurat.rda")
cells


#list of marker genes in the study 
marker_genes <- read.csv("/projectnb/bf528/project_4_scrnaseq/GSM2230760_marker_genes.csv")

#feature selection, identify highly variable features
cells <- FindVariableFeatures(cells, selection.method = "vst", nfeatures = 2000)

#top 10 most variable genes
top10 <- head(VariableFeatures(cells), 10)
top10

#plot variable features 
plot1 <- VariableFeaturePlot(cells)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
plot1

#scale the data
all_genes <- rownames(cells)
cells <- ScaleData(cells, features = all_genes)

#linear dimensional reduction on the scaled data, variable features as input 
cells <- RunPCA(cells, features = VariableFeatures(object = cells))

#examine and visualize PCA results 
VizDimLoadings(cells, dims = 1:2, reduction = "pca")

#visualize clusters by cluster and genes 
png("pca.png")
DimPlot(cells, reduction = "pca")

png("features.png")
FeaturePlot(cells, features = all_genes)

#clustering 
cells <- FindNeighbors(cells, dims = 1:10)
cells <- FindClusters(cells, resolution = 0.5)

#dimensional reduction
cells <- RunUMAP(cells, dims = 1:10)
cells <- RunTSNE(cells, dims = 1:10)
DimPlot(cells, reduction = "umap")
#DimPlot(cells, reduction = "")


#visualizing with violin plots 
png("violin_plots.png")
VlnPlot(cells, features = c("GCG","INS", "PPY", "SST", "GHRL", "PRSS1", "KRT19", "SPARC", "VWF", "RGS5", "PDGFRA", "SOX10", "SDS", "TPSAB1", "TRAC"))

#raw counts 
png("violin_plots_raw_counts.png")
VlnPlot(cells, features = c("GCG", "INS", "PPY", "SST", "GHRL", "PRSS1", "KRT19", "SPARC", "VWF", "RGS5", "PDGFRA", "SOX10", "SDS", "TPSAB1", "TRAC"), slot = "counts", log = TRUE)

#feature plot 
png("Featureplot.png")
FeaturePlot(cells, features = c("INS"))

#identify the marker genes for each cluster
#find all markers for every cluster compared to all remaining cells, 
cells_markers <- FindAllMarkers(cells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_csv("cells_markers", "cells_markers.csv")
#different method of clustering cell types 
top10_gene <- cells_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_gene
write.csv(top10_gene, "/projectnb/bf528/users/group_3/Project_4/analyst/top10_gene.csv")

#top 5 genes 
top5_gene <- cells_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5_gene


#violin plots of cells in the study
#alpha 
VlnPlot(cells, features = c("GCG")) #consistently expressed throughout
#beta and acinar 
VlnPlot(cells, features = c("INS"))
#Delta
VlnPlot(cells, features = c("SST"))
#gamma- do not have, not highly expressed in any of the clusters
VlnPlot(cells, features = c("PPY"))
#epsilon 
VlnPlot(cells, features = c("GHRL"))
#acinar 
VlnPlot(cells, features = c("CPA1"))
#ductal
VlnPlot(cells, features = c("KRT19"))
#quiescent stellate 
VlnPlot(cells, features = c("RGS5"))
#activated stellate 
VlnPlot(cells, features = c("PDGFRA"))
#endothelial cell
VlnPlot(cells, features = c("VWF"))
#macrophage 
VlnPlot(cells, features = c("SDS"))
#mast cell 
VlnPlot(cells, features = c("TPSAB1"))
#cytotoxic
VlnPlot(cells, features = c("TRCA"))


#top 5 marker genes in each cluster

#cluster 0- delta 
cluster0.markers <- FindMarkers(cells, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 10)

#cluster 1- beta 
cluster1.markers <- FindMarkers(cells, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 10)

#cluster 2 
cluster2.markers <- FindMarkers(cells, ident.1 = 2, min.pct = 0.25)
head(cluster2.markers, n = 10)

#cluster 3 
cluster3.markers <- FindMarkers(cells, ident.1 = 3, min.pct = 0.25)
head(cluster3.markers, n = 10)

#cluster 4
cluster4.markers <- FindMarkers(cells, ident.1 = 4, min.pct = 0.25)
head(cluster4.markers, n = 10)

#Cluster 5
cluster5.markers <- FindMarkers(cells, ident.1 = 5, min.pct = 0.25)
head(cluster5.markers, n = 10)

#Cluster 6
cluster6.markers <- FindMarkers(cells, ident.1 = 6, min.pct = 0.25)
head(cluster6.markers, n = 10)

#Cluster 7 
cluster7.markers <- FindMarkers(cells, ident.1 = 7, min.pct = 0.25)
head(cluster7.markers, n = 10)

#Cluster 8
cluster8.markers <- FindMarkers(cells, ident.1 = 8, min.pct = 0.25)
head(cluster8.markers, n = 10)

#Cluster 9
cluster9.markers <- FindMarkers(cells, ident.1 = 9, min.pct = 0.25)
head(cluster9.markers, n = 10)

#cluster 10
cluster10.markers <- FindMarkers(cells, ident.1 = 10, min.pct = 0.25)
head(cluster10.markers, n = 10)

#cluster 11
cluster11.markers <- FindMarkers(cells, ident.1 = 11, min.pct = 0.25)
head(cluster11.markers, n = 10)

#Cluster 12-
cluster12.markers <- FindMarkers(cells, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 10)

#Feature plots 
#Genes highly expressed

FeaturePlot(cells, features = c("SST"))
FeaturePlot(cells, features = c("GCG"))
FeaturePlot(cells, features = c("INS"))
FeaturePlot(cells, features = c("PPY"))

#Genes not expressed 
FeaturePlot(cells, features = c("TRAC"))
FeaturePlot(cells, features = c("TPSAB1"))
FeaturePlot(cells, features = c("PDGFRA"))
FeaturePlot(cells, features = c("VWF"))



#find markers distinguishing 1 from 8
cluster18.markers <- FindMarkers(cells, ident.1 = 1, ident.2 = 8, min.pct = 0.25)
head(cluster18.markers, n = 5)

#find markers distinguishing 3 from 4
cluster34.markers <- FindMarkers(cells, ident.1 = 3, ident.2 = 4, min.pct = 0.25)
head(cluster34.markers, n = 5)

#find markers distinguishing 0 from 2
cluster02.markers <- FindMarkers(cells, ident.1 = 0, ident.2 = 2, min.pct = 0.25)
head(cluster02.markers, n = 5)

#find markers distinguishing 4 from 6
cluster46.markers <- FindMarkers(cells, ident.1 = 4, ident.2 = 6, min.pct = 0.25)
head(cluster46.markers, n = 5)

#find markers distinguishing 6 from 7
cluster67.markers <- FindMarkers(cells, ident.1 = 6, ident.2 = 7, min.pct = 0.25)
head(cluster67.markers, n = 5)

# assign to cell types based on gene clusters 

#visualize cell types according to cluster 
#cell types mentioned in the paper: alpha (GCG), beta (INS), delta (SST), 
#Gamma (PPY), epsilon (GHRL), acinar (CPA1), ductal (KRT19), quiescent stellate (RGS5), 
#activated stellate(PDGFRA), endothelial (VWF), macrophage (SDS), mast(TPSAB1), cytoxic T (TRAC), Schwann (SOX10)

#assign cell types to clusters based on marker genes

new_cluster_ids <- c("Delta", "Beta", "Acinar", "Stellate", "Alpha", "Ductal", "Gamma", "Mesenchymal", "Beta", "Ductal", "Mesenchymal", "Macrophage", "Ductal")
names(new_cluster_ids) <- levels(cells)
cells <- RenameIdents(cells, new_cluster_ids)

#plot of clustered cells
png("dim_plot.png")
DimPlot(cells, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#visualize the top marker genes in the study 

png("Clustered_heatmap.png")
DoHeatmap(cells, features = top5_gene$gene) + NoLegend()

#additional heatmap
DoHeatmap(cells, features = top10_gene$gene) + NoLegend()

#saving seurat object 
saveRDS(cells, file = "analyst_results.rds")

#CSV of marker genes 
write_csv("marker_genes", "marker_genes.csv")
