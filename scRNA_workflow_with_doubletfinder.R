# This script is written to implement scRNA seq analysis 
# The workflow includes merge all samples for initial QC
# and subset for implementing DoubletFinders and merged again.

# Set working directory
# This script and data dir sit in the working directory
setwd("/data/regmib2/bio_data_mining/papa_sc/sc_workflow")

# load R package
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)

# Set seed for reproducibility
set.seed(1234)

# Organize data and create SeuratObject
dirs <- list.dirs(path="data", recursive = F, full.names = F)
for (x in dirs){
    name <- gsub('_filtered_feature_bc_matrix', '', x)
    cts <- ReadMtx(mtx = paste0('data/', x, '/matrix.mtx.gz'),
    features = paste0('data/', x, '/features.tsv.gz'),
    cells = paste0('data/', x, '/barcodes.tsv.gz'))
    assign(name, CreateSeuratObject(counts = cts))
}

# Merge all Seurat objects
merge_seurat <- merge(GM001_PRE, y = c(GM002_PRE, GM003_POST, GM004_POST),
    add.cell.ids = c("GM001_PRE", "GM002_PRE", "GM003_POST", "GM004_POST"),
    project = 'PAPA')
merge_seurat$sample <- rownames(merge_seurat@meta.data)
merge_seurat@meta.data <- separate(merge_seurat@meta.data, col = 'sample', 
    into = c("Patient","Type", "Barcode"), sep = "_")

# Calculate mitochondrial percentage
merge_seurat$percent_mito <- PercentageFeatureSet(merge_seurat, pattern = '^MT-')
merge_seurat$percent_ribo <- PercentageFeatureSet(merge_seurat, pattern = '^RP[SL]')
merge_seurat$percent_hb <- PercentageFeatureSet(merge_seurat, pattern = '^HB[^(P)]')
merge_seurat$percent_plat <- PercentageFeatureSet(merge_seurat, pattern = 'PECAM1|PF4')

# Plot QC
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(merge_seurat, group.by = "Patient", features = feats, pt.size = 0.1, ncol = 2) +
    NoLegend()

plot1 <- FeatureScatter(merge_seurat, "nCount_RNA", "percent_mito", group.by = "Patient", pt.size = 0.5)
plot2 <- FeatureScatter(merge_seurat, "nCount_RNA", "nFeature_RNA", group.by = "Patient", pt.size = 0.5)
plot1 + plot2

# Detection-based Filtering
selected_c <- WhichCells(merge_seurat, expression = nFeature_RNA > 200)
selected_f <- rownames(merge_seurat)[Matrix::rowSums(merge_seurat) > 3]

data.filt <- subset(merge_seurat, features = selected_f, cells = selected_c)
dim(data.filt)

# Check genes contribute the most to reads. # This might not be working
par(mar = c(2, 8, 2, 1))
C <- data.filt@assays$RNA@counts
C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
most_expressed_C <- order(apply(C, 1, median), decreasing = T)[20:1]

# boxplot the percentage of counts per gene.
boxplot(as.matrix(Matrix::t(C[most_expressed_C, ])), cex = 0.1, las = 1, 
    xlab = "% total count per cell", col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)

# Mito/Ribo filtering at cell level, and subset the object to only keep those cells
selected_mito <- WhichCells(data.filt, expression = percent_mito < 15)
selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.05)
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)

# Plot QC after filtering
feats_a <- c("nFeature_RNA", "nCount_RNA", "percent_mito", 
    "percent_ribo", "percent_hb", "percent_plat")
VlnPlot(data.filt, group.by = "Patient", features = feats_a, pt.size = 0.1, ncol = 2) +
    NoLegend()

plot3 <- FeatureScatter(data.filt, "nCount_RNA", "percent_mito", group.by = "Patient",
    pt.size = 0.5)
plot4 <- FeatureScatter(data.filt, "nCount_RNA", "nFeature_RNA", group.by = "Patient", 
    pt.size = 0.5)
plot3 + plot4

# Filter Mitocondrial & ribosomal genes
# high detection of MALAT1 may be a technical issue
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
data.filt <- data.filt[ ! grepl('^RP[SL]', rownames(data.filt)), ]

# Filter Hemoglobin gene (optional if that is a problem on the data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]

# check sample gender
library(biomaRt) # ask for the file
genes.table = read.csv("/data/jiangk/projects/BMDS/external/biomaRt/genes.file")
genes.table <- genes.table[genes.table$external_gene_name %in% rownames(data.filt),]
chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]
data.filt$pct_chrY = 
    colSums(data.filt@assays$RNA@counts[chrY.gene, ])/colSums(data.filt@assays$RNA@counts)

# Plot XIST expression vs chrY proportion
FeatureScatter(data.filt, feature1 = "XIST", feature2 = "pct_chrY")
VlnPlot(data.filt, features = c("XIST", "pct_chrY"))

# Calculate cell-cycle scores
data.filt <- NormalizeData(data.filt)
data.filt <- CellCycleScoring(object = data.filt,
    g2m.features = cc.genes$g2m.genes,
    s.features = cc.genes$s.genes)

# Warning: The following features are not present in the object: MLF1IP, not searching for symbol 
# synonyms
# Warning: The following features are not present in the object: FAM64A, HN1, not searching for 
# symbol synonyms
VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "Patient",
    ncol = 2, pt.size = 0.1)

# Subset samples to implement DoubletFinder

require(DoubletFinder)

obj.list <- SplitObject(data.filt, split.by = "Patient")
for(i in 1:length(obj.list))
    {print(paste0("sample", i))
    object.sample <- NormalizeData(obj.list[[i]])
    object.sample <- FindVariableFeatures(object.sample, verbose = F)
    object.sample <- ScaleData(object.sample)
    object.sample <- RunPCA(object.sample, verbose = F, unpcs = 20)
    object.sample <- RunUMAP(object.sample, dims = 1:20, verboose = F)
  
    sweep.list <- paramSweep_v3(object.sample, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.list)
    bcmvn <- find.pK(sweep.stats)
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    nExp <- round(ncol(object.sample) * 0.13)
    list.Doublet.Finder.Parameter[i] <- paste0("sample  ", i,"  ", optimal.pk)
# Implement DoubletFinder with optimal pK paramete
# Note pN shoould be users defined
    object.sample <- doubletFinder_v3(seu = object.sample, pN = 0.25, pK = optimal.pk, 
    nExp = nExp, PCs = 1:20)
  
# For each sample, DoubletFinder adds column with different column name, therefore column 
# names should be changed
    metadata <- object.sample@meta.data
    colnames(metadata)[ncol(metadata)] <- "doublet_type"
    colnames(metadata)[(ncol(metadata)-1)] <- "doublet_estimate"
    object.sample@meta.data <- metadata
    obj.list[[i]] <- object.sample
  }

# Combine the seurat object list
combined.filt.doublets <- obj.list[[1]]
for (i in 2:length(obj.list)){
    combined.filt.doublets <- merge(x = combined.filt.doublets, y=obj.list[[i]])
}

# Visualize doublets
combined.filt.doublets <- FindVariableFeatures(combined.filt.doublets, verbose = F)
combined.filt.doublets <- ScaleData(combined.filt.doublets)
combined.filt.doublets <- RunPCA(combined.filt.doublets, verbose = F, npcs = 20)
combined.filt.doublets <- RunUMAP(combined.filt.doublets, dims = 1:20, verbose = F)
DimPlot(combined.filt.doublets, reduction = 'umap', group.by = "doublet_type")
pre_doub_filt_dimension <- dim(combined.filt.doublets)

# Filter out doublets
DF.name <- colnames(combined.filt.doublets@meta.data)[grepl("doublet_type", 
    colnames(combined.filt.doublets@meta.data))]
combined.doublet.free <- combined.filt.doublets[, combined.filt.doublets@meta.data[, DF.name] == "Singlet"]
post.doub.filt.dimension <- dim(combined.doublet.free)

# Annotation
library(monocle)
library(SeuratWrappers)
library(harmony)
library(patchwork)
library(scran)
library(cowplot)

set.seed(348)
# Run the standard workflow for visualization and clustering
combined <- NormalizeData(object = combined.doublet.free)
combined <- FindVariableFeatures(object = combined)
combined <- ScaleData(object = combined)
combined <- RunPCA(object = combined)
ElbowPlot(combined)
combined <- FindNeighbors(object = combined, dims = 1:20)
combined <- FindClusters(object = combined)
combined <- RunUMAP(object = combined, dims = 1:20)

DefaultAssay(combined) <- "RNA"
Combined.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, 
    logfc.threshold = 0.25)
top15 <- Combined.markers %>% group_by(cluster) %>% top_n(n = 15, wt = avg_log2FC)

# Plot of genes of interest
myfeatures <- c("CD3E", "CD4", "CD8A", "NKG7", "GNLY", "MS4A1", "CD14", "LYZ", "MS4A7",
                "FCGR3A", "CST3", "FCER1A")
plot_list <- list()
for (i in myfeatures) {
    plot_list[[i]] <- FeaturePlot(combined, reduction = "umap", dims = 1:2,
    features = i, split.by = "Patient", min.cutoff = 'q10',
    ncol = 2,label = FALSE,
    order = T) + NoLegend() + NoAxes() + NoGrid()
}

plot_grid(ncol = 2, plotlist = plot_list)


# Annotation using SingleR
# created singlel cell experiment
combined.sce <- as.SingleCellExperiment(combined, assay = "RNA")

library(SingleR)
library(celldex)
library(pheatmap)

# now label using singleR and celldex (requires an internet connection to connect to ExperimentHub)
#version below: 2022-10-31
HPCA.data <- HumanPrimaryCellAtlasData(ensembl = FALSE)
DICE.data <- DatabaseImmuneCellExpressionData(ensembl = FALSE)
Monaco.data <- MonacoImmuneData(ensembl = FALSE)

# Combined.sce <- as.SingleCellExperiment(Combined, assay = "RNA")
predictions.dice <- SingleR(test=combined.sce, assay.type.test=1,
    ref=DICE.data, labels=DICE.data$label.fine)

predictions.monaco <- SingleR(test=combined.sce, assay.type.test=1,
    ref=Monaco.data, labels=Monaco.data$label.fine)

predictions.hpca <- SingleR(test=combined.sce, assay.type.test=1,
    ref=HPCA.data, labels=HPCA.data$label.fine)

# Add back to singleCellExperiment object (or Seurat objects)
combined.sce[["SingleR.labels.dice"]] <- predictions.dice$labels
combined.sce[["SingleR.labels.monaco"]] <- predictions.monaco$labels
combined.sce[["SingleR.labels.hpca"]] <- predictions.hpca$labels

# Cell type to Seurat object # there could be difference between .labels and pruned.labels, 
# I just put purned labels here
combined$celltype_dice <- predictions.dice$labels
combined$celltype_monaco <- predictions.monaco$labels
combined$celltype_hpca <- predictions.hpca$labels

# cell type to seurat object
combined$celltype_dice_purned_labels <- predictions.dice$pruned.labels
combined$celltype_monaco_pruned_labels <- predictions.monaco$pruned.labels
combined$celltype_hpca_purned_labels <- predictions.hpca$pruned.labels

combined <- RunUMAP(object = combined, dims = 1:20)

# visualization
# Examine and visualize PCA results a few different ways
DimPlot(combined, reduction = "umap", split.by = "Patient", label = TRUE)

# 
combined2 <- as.Seurat(combined.sce, counts = NULL)

DimPlot(combined2, reduction = "UMAP",
    group.by = "SingleR.labels.dice",
    label = TRUE,label.size = 3)

DimPlot(combined2, reduction = "UMAP",
    group.by = "SingleR.labels.monaco",
    label = TRUE,label.size = 3)

DimPlot(combined2, reduction = "UMAP",
    group.by = "SingleR.labels.hpca",
    label = TRUE,label.size = 3)

# Assign cell type
new.cluster.ids <- c("CD14+CD16-Mono", "CD8 Naive", "NK_CD56+", "CD4 Naive", "NK_CD56-", "Treg", 
    "?CD4 effector","CD14+CD16+Mono","?NA","DC_CD74+","B cells","?NA_T", "Platelets", 
    "Platelets", "Ery cells", "DC_CD74-", "Mono/DC_IL1B+", "NK3", "NK4", "Ery cells", "NA")

names(new.cluster.ids) <- levels(combined2)
combined2 <- RenameIdents(combined2, new.cluster.ids)
DimPlot(combined2, reduction = "UMAP",
    split.by = "Patient",
    label = TRUE)

save.image(file = "IMAGE_scRNA_with_doubletfinder.RData")













