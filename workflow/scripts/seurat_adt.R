library(Seurat)
library(dplyr)
library(MAST)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
# library(scater)


## ----Input, include=FALSE------------------------------------------------------------
args<-commandArgs(TRUE);
workdir <- args[1]
data_path <- args[2]
sample <- args[3]
genome <- args[4]


## ----Load Data-----------------------------------------------------------------------
rdata <- Read10X(data_path)


## ----Create Working Directory, include=FALSE-----------------------------------------
dir.create(workdir, recursive=TRUE)


## ----Change Working Directory, include=FALSE-----------------------------------------
setwd(workdir)


## ----Create Seurat Object------------------------------------------------------------
seur <- CreateSeuratObject(counts=rdata$`Gene Expression`, project = sample)

adt_assay <- CreateAssayObject(counts=rdata$`Antibody Capture`[grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE, invert=TRUE),])
seur[["ADT"]] <- adt_assay

hashtag = FALSE
if (length(grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE)) > 0) {
  hto_assay <- CreateAssayObject(counts=rdata$`Antibody Capture`[grep('hashtag', rownames(rdata$`Antibody Capture`), value=TRUE, ignore.case=TRUE),])
  seur[["HTO"]] <- hto_assay
  hashtag = TRUE
}


## ----Calculate Mitochondrial Percentage----------------------------------------------
if (genome == "mm10") {
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^mt-")
} else if (genome == "hg19" | genome == "hg38") {
  seur[["percent.mito"]] <- PercentageFeatureSet(seur, pattern="^MT-")
}


## ----PreFilter Gene Plot Plot, echo=FALSE, message=FALSE, results="hide", fig.width=10, fig.height=5----
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito") + NoLegend()
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + NoLegend()

png("PreFilter_Gene_Plot.png", height=5, width=10, units='in', res=300)
plot1+plot2
dev.off()


## ----Get Filtering Thresholds--------------------------------------------------------
#Filtering
thresh <- list()
thresh['nFeature_RNA_low'] <- median(log1p(seur$nFeature_RNA)) - 3*mad(log1p(seur$nFeature_RNA))
thresh['nFeature_RNA_low'] <- median(log1p(seur$nFeature_RNA)) - 3*mad(log1p(seur$nFeature_RNA))
thresh['nFeature_RNA_high'] <- median(log1p(seur$nFeature_RNA)) + 3*mad(log1p(seur$nFeature_RNA))
thresh['nCount_RNA_low'] <- median(log1p(seur$nCount_RNA)) - 3*mad(log1p(seur$nCount_RNA))
thresh['nCount_RNA_high'] <- median(log1p(seur$nCount_RNA)) + 3*mad(log1p(seur$nCount_RNA))
thresh['mt_high'] = median(log1p(seur$percent.mito)) + 3*mad(log1p(seur$percent.mito))

thresh['nFeature_ADT_low'] <- median(log1p(seur$nFeature_ADT)) - 3*mad(log1p(seur$nFeature_ADT))
thresh['nFeature_ADT_high'] <- median(log1p(seur$nFeature_ADT)) + 3*mad(log1p(seur$nFeature_ADT))
thresh['nCount_ADT_low'] <- median(log1p(seur$nCount_ADT)) - 3*mad(log1p(seur$nCount_ADT))
thresh['nCount_ADT_high'] <- median(log1p(seur$nCount_ADT)) + 3*mad(log1p(seur$nCount_ADT))

if (hashtag) {
  thresh['nCount_HTO_high'] <- median(log1p(seur$nCount_HTO)) + 3*mad(log1p(seur$nCount_HTO))
}


## ----Get Cells to Remove-------------------------------------------------------------
cellsToRemove.Feature_RNA <- colnames(seur)[which(log1p(seur$nFeature_RNA) < thresh['nFeature_RNA_low'] | log1p(seur$nFeature_RNA) > thresh['nFeature_RNA_high'])]
cellsToRemove.Count_RNA <- colnames(seur)[which(log1p(seur$nCount_RNA) < thresh['nCount_RNA_low'] | log1p(seur$nCount_RNA) > thresh['nCount_RNA_high'])]
cellsToRemove.mito <- colnames(seur)[which(log1p(seur$percent.mito) > thresh['mt_high'])]

#cellsToRemove.Feature_ADT <- colnames(seur)[which(log1p(seur$nFeature_ADT) < nFeature_ADT_low | log1p(seur$nFeature_ADT) > nFeature_ADT_high)]
cellsToRemove.Count_ADT <- colnames(seur)[which(log1p(seur$nCount_ADT) > thresh['nCount_ADT_high'])]

cellsToRemove.Count_HTO <- c()
if (hashtag) {
  cellsToRemove.Count_HTO <- colnames(seur)[which(log1p(seur$nCount_HTO) > thresh['nCount_HTO_high'])]
}

thresh['numCellsRemove'] <- length(unique(c(cellsToRemove.Feature_RNA, cellsToRemove.Count_RNA, cellsToRemove.mito, cellsToRemove.Count_ADT, cellsToRemove.Count_HTO)))

write.table(thresh, 'cell_filter_info.csv', quote=FALSE, row.names=FALSE, sep=',')



## ----Pre-Filter RNA Violin Plot, warning=FALSE, message=FALSE, results="hide", fig.height=7----
##Pre-Filter Violin Plot
plot1 <- VlnPlot(seur, features = c("nFeature_RNA")) + NoLegend() + geom_hline(yintercept=exp(thresh$nFeature_RNA_low)-1,linetype="dashed") + geom_hline(yintercept=exp(thresh$nFeature_RNA_high)-1,linetype="dashed")
plot2 <- VlnPlot(seur, features = "nCount_RNA") + NoLegend() + geom_hline(yintercept=exp(thresh$nCount_RNA_low)-1,linetype="dashed") + geom_hline(yintercept=exp(thresh$nCount_RNA_high)-1,linetype="dashed")
plot3 <- VlnPlot(seur, features = "percent.mito", ncol = 3) + NoLegend() + geom_hline(yintercept=exp(thresh$mt_high)-1,linetype="dashed")

png("PreFilter_VlnPlot_RNA.png", height=7, width=7, units='in', res=300)
grid.arrange(plot1, plot2, plot3, nrow=1)
dev.off()


## ----Pre-Filter ADT Violin Plot, warning=FALSE, message=FALSE, results="hide", fig.width=5, fig.height=7, out.width="70%"----
plot1 <- VlnPlot(seur, features = c("nFeature_ADT")) + NoLegend() 
plot2 <- VlnPlot(seur, features = "nCount_ADT") + NoLegend() + geom_hline(yintercept=exp(thresh$nCount_ADT_high)-1,linetype="dashed") 

png("PreFilter_VlnPlot_ADT.png", height=7, width=5, units='in', res=300)
grid.arrange(plot1, plot2, nrow=1)
dev.off()


## ----Pre-Filter HTO Violin Plot, warning=FALSE, message=FALSE, results='hide', fig.width=5, fig.height=7, out.width="70%", eval=hashtag----

if (hashtag) {
  plot1 <- VlnPlot(seur, features = c("nFeature_HTO")) + NoLegend() 
  plot2 <- VlnPlot(seur, features = "nCount_HTO") + NoLegend() + geom_hline(yintercept=exp(thresh$nCount_HTO_high)-1,linetype="dashed") 

  png("PreFilter_VlnPlot_HTO.png", height=7, width=5, units='in', res=300)
  grid.arrange(plot1, plot2, nrow=1)
  dev.off()
}


## ----Subset Cells--------------------------------------------------------------------
seur <- subset(seur, cells = unique(c(cellsToRemove.Feature_RNA, cellsToRemove.Count_RNA, cellsToRemove.mito, cellsToRemove.Count_ADT, cellsToRemove.Count_HTO)), inver=T)


## ----Demuxlet---------------

demuxbest <- read.table('demuxlet/output/S2/S2.best', sep='\t', header=TRUE)

rownames(demuxbest) <- demuxbest$BARCODE

seur <- AddMetaData(seur, metadata = demuxbest[colnames(seur),])

seur <- AddMetaData(seur, metadata = sapply(seur$BEST, function(x) strsplit(x, '-')[[1]][[1]]), col.name='DROPLET.TYPE')

write.table(table(seur$DROPLET.TYPE), 'demuxlet_droplet_counts.csv', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')

seur.full <- seur
seur.AMB <- subset(seur, subset = DROPLET.TYPE == "AMB")
seur <- subset(seur, subset = DROPLET.TYPE == "SNG")

write.table(table(seur$BEST), 'demuxlet_singlet_counts.csv', row.names=FALSE, col.names=FALSE, quote=FALSE, sep=',')


## ----Post-Filter Gene Plot, echo=FALSE, warning=FALSE, message=FALSE, results="hide", fig.width=10, fig.height=5----
#Post-Filter Plots
plot1 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "percent.mito")
plot2 <- FeatureScatter(seur, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
png("PostFilter_Gene_Plot.png", height=5, width=10, units='in', res=300)
plot1+plot2
dev.off()


## ----Post-Filter RNA Violin Plot, results="hide", fig.height=7-----------------------
png("PostFilter_VlnPlot_RNA.png", height=7, width=7, units='in', res=300)
VlnPlot(seur, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
dev.off()


## ----Post-Filter ADT Violin Plot, results="hide", fig.width=5, fig.height=7, out.width="70%"----
png("PostFilter_VlnPlot_ADT.png", height=7, width=5, units='in', res=300)
VlnPlot(seur, features = c("nFeature_ADT", "nCount_ADT"), ncol = 2)
dev.off()


## ----Post-Filter HTO Violin Plot, warning=FALSE, message=FALSE, results='hide', fig.width=5, fig.height=7, out.width="70%", eval=hashtag----

if (hashtag) {
  png("PostFilter_VlnPlot_HTO.png", height=7, width=5, units='in', res=300)
  print(VlnPlot(seur, features = c("nFeature_HTO", "nCount_HTO"), ncol = 2))
  dev.off()
}


## ----Normalize RNA Data, message=FALSE-----------------------------------------------
#Normalize RNA data
seur <- NormalizeData(seur) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
seur <- FindNeighbors(seur, dims = 1:30)
seur <- FindClusters(seur, resolution = 0.8, algorithm=3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                      reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')


## ----Normalize HTO Data, warning=FALSE, message=FALSE, eval=hashtag------------------
if (hashtag) {
  seur <- NormalizeData(seur, assay='HTO', normalization.method = 'CLR', margin = 2)
  seur <- HTODemux(seur, assay = "HTO", positive.quantile = 0.99)
  seur$hash.ID <- factor(seur$hash.ID, levels=levels(seur$hash.ID)[order(levels(seur$hash.ID))])

## ----HTO Ridge Plots, results='show', eval=hashtag, fig.width=12, fig.height=9-------
  png('HTO_Ridge_Plot.png', units='in', width=12, height=9, res=300)
  for (i in seq(1,length(rownames(seur[["HTO"]])), by=25)) {
    print(RidgePlot(seur, sort(rownames(seur[['HTO']]))[i:min(i+24,length(rownames(seur[['HTO']])))], assay="HTO", ncol=min(5, ceiling(sqrt(length(rownames(seur[['HTO']]))-(i-1)))), group.by='hash.ID'))
  }
}


## ----RNA UMAP Plot, echo=FALSE, warning=FALSE, message=FALSE, results="hide"---------
png("UMAP_RNA.png", width=1800, height=1600, res = 300)
DimPlot(seur, reduction='rna.umap', group.by='RNA_snn_res.0.8', label = TRUE) + ggtitle("RNA")
dev.off()


## ----Normalize ADT Data, message=FALSE-----------------------------------------------
#Normalize ADT data
DefaultAssay(seur) <- 'ADT'
# we will use all ADT features for dimensional reduction
# we set a dimensional reduction name to avoid overwriting the 
VariableFeatures(seur) <- rownames(seur[["ADT"]])
seur <- NormalizeData(seur, normalization.method = 'CLR', margin = 2) %>% 
  ScaleData() %>% RunPCA(reduction.name = 'apca')
seur <- FindNeighbors(seur, dims = 1:min(length(rownames(seur[['ADT']]))-1, 20), reduction = "apca")
seur <- FindClusters(seur, graph.name = "ADT_snn", algorithm = 3, verbose = FALSE)
seur <- RunUMAP(seur, reduction = 'apca', dims = 1:min(length(rownames(seur[['ADT']]))-1, 20), assay = 'ADT', 
                reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')


## ----ADT UMAP Plot, echo=FALSE, warning=FALSE, message=FALSE, results="hide"---------
png("UMAP_ADT.png", width=1800, height=1600, res = 300)
DimPlot(seur, reduction='adt.umap', label = TRUE) + ggtitle("RNA")
dev.off()


## ----Weighted Nearest Neighbors, message=FALSE---------------------------------------
#Get nearest neighbors based on both RNA and ADT data
seur <- FindMultiModalNeighbors(
  seur, reduction.list = list("pca", "apca"), 
  dims.list = list(1:30, 1:dim(seur@reductions$apca)[[2]]), modality.weight.name = "RNA.weight"
)


#UMAP based on the multi-modal weighted graph
seur <- RunUMAP(seur, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
seur <- FindClusters(seur, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


## ----Multi-Modal UMAP Plot, echo=FALSE, warning=FALSE, message=FALSE, result='hide'----
p1 <- DimPlot(seur, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) #+ NoLegend()

png("UMAP_MultiModal_WNN.png", width=1800, height=1600, res = 300)
p1 + ggtitle("Multi-Modal Weighted Graph")
dev.off()

## ----Multi-Modal Modality Weights, echo=FALSE, warning=FALSE, message=FALSE, results='hide', fig.width=10, fig.height=7----
png("MultiModal_WNN_Weights.png", width=3600, height=2400, res = 300)
(VlnPlot(seur, features = "RNA.weight", group.by = 'wsnn_res.0.8', pt.size = 0.1) + NoLegend()) / (VlnPlot(seur, features = "ADT.weight", group.by = 'wsnn_res.0.8', pt.size = 0.1) + NoLegend())
dev.off()

## ----RNA and ADT UMAP Plot, echo=FALSE, warning=FALSE, message=FALSE, result='hide', fig.width=12----
#UMAP based on RNA and ADT data individually

p3 <- DimPlot(seur, reduction = 'rna.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5) + NoLegend() + ggtitle("RNA")
p4 <- DimPlot(seur, reduction = 'adt.umap', label = TRUE, 
              repel = TRUE, label.size = 2.5)  + ggtitle("ADT") #+ NoLegend()
png("UMAP_MultiModal_RNA_ADT.png", width=3600, height=1600, res = 300)
p3 + p4
dev.off()


## ----RidgePlots of ADT, echo=FALSE, warning=FALSE, message=FALSE, result='hide'------
adt_thresh = 0
suppressWarnings(dir.create('ADT'))
#pdf('./Ridgeplots.pdf', width=21, height=21)
count <- 1
for (i in seq(1,length(names(which(rowSums(seur[['ADT']]) > adt_thresh))), by=25)) {
  png(paste0("ADT/ADT_Ridgeplots_", count, ".png"), width=21, height=21, res=300, units='in')
  print(RidgePlot(seur, sort(names(which(rowSums(seur[['ADT']]) > adt_thresh)))[i:min(i+24,length(rownames(seur[['ADT']])))], assay="ADT", ncol=5))
  dev.off()
  count <- count + 1
}
#dev.off()

write.table(adt_thresh, 'ADT_threshold.txt', col.names = FALSE, row.names=FALSE)

if (length(which(rowSums(seur[['ADT']]) <= adt_thresh)) > 0){
  write.table(sort(names(which(rowSums(seur[['ADT']]) <= adt_thresh))), 'ADT_excluded.csv', quote=FALSE, sep=',', col.names = FALSE, row.names = FALSE)
}

## ----Prepare for Heatmap, include=FALSE----------------------------------------------
require(scales)
identities <- levels(seur$wsnn_res.0.8)
DefaultAssay(seur) <- 'RNA'

my_color_palette <- list(hue_pal()(length(identities)))

cols <- as.character(c((1:length(my_color_palette[[1]]))) - 1)
colors <- list(seurat_clusters=array(my_color_palette[[1]], dim = length(my_color_palette[[1]]), dimnames = list(cols)))

pbmc_subset.sce <- as.SingleCellExperiment(seur)
suppressWarnings(dir.create('DEG'))


## ----DEG Calculation, message=FALSE, results='asis'----------------------------------
cat("## Differential Expression - MAST {.tabset}")


# 1 cluster vs all other clusters
 for (cluster in c(0:(length(identities) - 1))) {
   try({
     Markers <- FindMarkers(seur,ident.1=cluster,test.use="MAST",only.pos = T)
     if (length(row.names(Markers) != 0)){
       Markers$sub <- Markers$pct.1-Markers$pct.2
       Markershead <- head(Markers, n=20)
       
       filename <- paste0('DEG/', sample, "_", cluster, "_markers.txt")
       write.table(cbind(Genes = rownames(Markers), Markers), filename, sep="\t", quote=FALSE, row.names=FALSE)
       
       Markers12 <- head(rownames(Markershead), n=12)
       filename <- paste0('DEG/', sample, "_", cluster,"_Features.png", sep="")
       png(filename , width=3200, height=3200, res = 300)
       print(FeaturePlot(seur, features = Markers12, ncol=3))
       dev.off()
     }
   })
 }


## ----Save RDS Object-----------------------------------------------------------------
saveRDS(seur, file="seur_cite_cluster.rds")


## ----Output Session Info, echo=TRUE, results='show'----------------------------------
sessionInfo()

