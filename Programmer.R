#'Packages needed:
#'(1) R server 4.2.3  (2)BiocManager::install(version = "3.16", "BioManager")
#'(3) BiocManager::install("tximeta") (4) BiocManager::install("org.Hs.eg.db") 
#'and then loading multiple library

library(patchwork)
library(tidyverse)
#' Alevin is a tool for estimating gene abundances in droplet-based single-cell dscRNA-seq data. 
#' Alevin-Seurat Connection
library(Seurat)
# for mapping Ensembl gene id to gene symbol (homo sapien)
library(reticulate)
library(umap)
library(tximeta)
library(SingleCellExperiment)
library(org.Hs.eg.db)

##' Preprocess the files in the projectnb toy data to be used in Seurat
#' path to the output directory of Alevin run of data
#files <- file.path("/projectnb/bf528/project_4_scrnaseq/GSM2230760__salmon_quant/alevin/quants_mat.gz")
files <- file.path("/projectnb/bf528/users/group_3/Project_4/Data_Curator/salmon_output/alevin/quants_mat.gz")
#file.exists(files) ---TRUE

#' Reading in the alevin quants quants
#txi <- tximport(files, type="alevin")
se <- tximeta(files, type="alevin", alevinArgs=list(filterBarcodes=TRUE))

#' tximeta returns a SummarizedExperiment. We can easily convert this object into a SingleCellExperiment,
#' which has specific slots designed for single-cell experiment data.
#' (SingleCellExperiment object is used widely across Bioconductor packages)
sce <- as(se, "SingleCellExperiment")

#'The data is now available as assays in sce
#'assayNames(sce)
#'## [1] "counts"

#'Then we can automatically add gene IDs, because tximeta knows the type of identifiers on the rows of the sce object:
#library(org.Hs.eg.db) # org pkg for Homo sapiens were used in this case
sce <- addIds(sce, "SYMBOL")
# mcols(sce) # could see the result

##'but first remove the ensembl id version num
ensembl_id <- rownames(sce)%>%
  as_vector()
ensembl_id <- gsub('\\.[0-9]*$', '', ensembl_id)
rownames(sce)<-ensembl_id

#Create a df to let the correct ensembl_id and corresponding symbol store in this df
id_df<-sce@rowRanges@elementMetadata$gene_id%>%as.data.frame()
colnames(id_df)[1] ="ensembl_id"
id_df["id"]<- gsub('\\.[0-9]*$', '', id_df[['ensembl_id']])
id_df['symbol']<-sce@rowRanges@elementMetadata$SYMBOL
id_filt<- id_df[!is.na(id_df$symbol),]
# Subset sce@assays@data$counts to keep only rows with id that match those in id_filt
rownames(sce@assays@data$counts)<-id_df$id
subset_sce <- sce@assays@data$counts[rownames(sce@assays@data$counts) %in% id_filt$id,]
# Create a named vector of rownames to replace
replacement <- setNames(id_filt$symbol, id_filt$id)
# Replace the rownames in subset_txi with their corresponding symbol values
rownames(subset_sce) <- replacement[rownames(subset_sce)]
# Combine the subset and remaining rows of scpc
sce@assays@data$counts<- rbind(subset_sce, sce@assays@data$counts[!rownames(sce@assays@data$counts) %in% id_filt$id,])
# Convert row names to vector and subset matrix
sce@assays@data$counts <- sce@assays@data$counts[!duplicated(rownames(sce@assays@data$counts)), ]




#'and we are good to go !! Cells after this has been taken from Seurat tutorial:
#'use count matric to create a Seurat Object (10x Genomics is a microfluidics-based method of single-cell RNA sequencing.)
cts <- sce@assays@data$counts

#Create Seurat Object
scpc <- CreateSeuratObject(cts, min.cells = 3, min.features = 200)

#' count matrix QC and selecting cells for further analysis
#' calculate mitochondrial QC metrics with the PercentageFeatureSet function
mt_genes <- rownames(sce)[as.logical(seqnames(sce) == "chrM")]
mt_genes <- mt_genes[mt_genes %in% rownames(GetAssayData(object = scpc, slot = "counts"))]
scpc[["percent.mt"]] <- PercentageFeatureSet(scpc, features= mt_genes)

# Set the resolution of the plot to 72 dpi
#png(filename = "my_plot.png", width = 900, height = 900, res = 250)

#' Visualize QC metrics as a violin plot to let us decide how we are gonna to filter the low-quality cells
# Generate the plot
VlnPlot(scpc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Save and close the plot
#dev.off()

##'Plot scatter plots to set the filter threshold
FeatureScatter(object = scpc, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = scpc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


# could use "View(scpc@meta.data)" to check all the values of "nFeature_RNA, nCount_RNA"
#' Since the majority of the dots fall in the lower region of the plot. 
#' and the violin plot is obscured by the dense concentration of dots, so I would like to have more specific numbers to assist me to decide the filter threshold
meta <- scpc@meta.data
nFeature_info <- summary(meta$nFeature_RNA) #Mean 656.5 #Median 289
nCount_info <- summary(meta$nCount_RNA) #Mean 1135 # Median 375
percent.mt_info <- summary(meta$percent.mt) #Mean 18.99 #Median 17.92
########But my mt_percentage is all zeros!!!!!

#turn the Feature_RNA column in scpc@meta.data to numeric data type
#scpc@meta.data$nFeature_RNA <- as.numeric(scpc@meta.data$nFeature_RNA)
scpc <- subset(scpc, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 &  percent.mt < 15)
#scpc <- subset(scpc, subset = nFeature_RNA > 200 & nFeature_RNA < 3500 &  percent.mt < 55)

#' Normalizing the data
scpc <- NormalizeData(scpc, normalization.method = "LogNormalize", scale.factor = 10000)

#identification of highly variable features (feature selection)
#scpc <- FindVariableFeatures(scpc, selection.method = "vst", nfeatures = 656)
scpc <- FindVariableFeatures(scpc, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scpc), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(scpc)
#LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0) #This line does not work!!

#'Scaling the data
all_genes <- rownames(scpc)
scpc <- ScaleData(scpc, features = all_genes)
#' Perform linear dimensional reduction
scpc <- RunPCA(scpc, features = VariableFeatures(object = scpc))
#'plot all of our projections
DimPlot(scpc, reduction = "pca")
#'elbow plot which simply quantitates the amount of variance captured in the different PCs
ElbowPlot(scpc) #there are fairly high amounts of information captured in the first 13 PCs. Taking somewhere between 12-17 PCs should therefore capture what we want to see.
#'easy exploration of the primary sources of heterogeneity in a dataset
DimHeatmap(scpc,dims=1:15, cells=500)


#before this you need to install the 

#' Cluster the cells
scpc <- FindNeighbors(scpc, dims = 1:10)
scpc <- FindClusters(scpc, resolution = 0.7)
#'Look at cluster IDs of the first 5 cells
head(Idents(scpc), 5)
cluster<-Idents(scpc)

#' Run non-linear dimensional reduction (UMAP/tSNE)
final_scpc <- RunUMAP(scpc, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(final_scpc, reduction = "umap",label = TRUE)

#' Count the number of cells in each cluster
cluster_tib <- cluster%>%as_tibble()
freq_table <- table(cluster_tib$value)%>%as.data.frame()
print(freq_table)
sum(freq_table['Freq'])

#Check cluster correctness
cluster1.markers <- FindMarkers(final_scpc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)

#Turn Seurat Object into RDS file
saveRDS(final_scpc, file = "/usr4/bf527/ivycwf/Documents/BF528/Project4/scpc_wo.rds")
saveRDS(final_scpc, file = "/projectnb/bf528/users/group_3/Project_4/Programmer/GSM2230760_seurat2.rds")







#Draw the stack barplot for illustrate the proportion of clusters in this sample
# Stacked + percent
freq_table <- freq_table %>% 
                mutate(sample = "human_sample")
# Define the number of colors you want
library(RColorBrewer)
nb.cols <- 14
mycolors <- colorRampPalette(brewer.pal(8, "Spectral"))(nb.cols)

# Sort the data frame by descending Freq values
freq_table <- freq_table[order(-freq_table$Freq),]

# Reorder Var1 factor levels based on sorted Freq values
ggplot(freq_table, aes(fill = reorder(Var1, Freq), y = Freq, x = sample)) +
  geom_bar(position = "fill", stat = "identity", width = 0.2, colour = "black", size = 0.3) +
  labs(fill = "Clusters") +
  scale_fill_manual(values = mycolors) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme_classic() +
  theme(axis.title.y = element_blank())






