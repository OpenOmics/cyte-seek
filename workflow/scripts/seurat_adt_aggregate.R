library(Seurat)
library(dplyr)
library(MAST)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)

library(future)
library(future.apply)
plan("multisession", workers = 8)
options(future.globals.maxSize = 30000 * 1024^2)

args<-commandArgs(TRUE);
workdir <- args[1]
genome <- args[2]

ref_samples <- args[3]
files <- args[4:length(args)]

data.list <- lapply(files, readRDS)

## ----Create Working Directory, include=FALSE-----------------------------------------
dir.create(workdir, recursive=TRUE)


## ----Change Working Directory, include=FALSE-----------------------------------------
setwd(workdir)

data.list <- lapply(X = data.list, FUN = function(x) {
  if ('SCT' %in% names(x)) {
    x[['SCT']] <- NULL
  }
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

#data.list <- lapply(X = data.list, FUN = SCTransform)

reference.features <- SelectIntegrationFeatures(object.list = data.list, nfeatures = 3000)

#data.list <- PrepSCTIntegration(object.list = data.list, anchor.features = reference.features)

data.list <- lapply(X = data.list, FUN = function(x) {
  x <- ScaleData(x, features = reference.features, verbose = FALSE)
  x <- RunPCA(x, features = reference.features, verbose = FALSE)
})

totalcells <- sum(sapply(data.list, function(x) dim(x)[[2]]))
future.seed=TRUE

ref_index <- sapply(strsplit(ref_samples, ',')[[1]], function(x) for (j in 1:length(data.list)) { if (x == data.list[[j]]@project.name) { return(j) }})
#reference.anchors = FindIntegrationAnchors(object.list = data.list, reduction = "rpca", dims = 1:30, anchor.features = reference.features, normalization.method = "SCT")
reference.anchors <- FindIntegrationAnchors(object.list = data.list, reference = ref_index, reduction = "rpca",
                                  dims = 1:50, normalization.method = "LogNormalize")

future.seed=TRUE
#combinedObj.integrated = IntegrateData(anchorset = reference.anchors, new.assay.name = "integratedRNA", dims = 1:30, normalization.method = "SCT")
combinedObj.integrated = IntegrateData(anchorset = reference.anchors, new.assay.name = "integratedRNA", dims = 1:50,  normalization.method = "LogNormalize")

for (i in grep('snn_res', colnames(combinedObj.integrated@meta.data), value=TRUE)) {
    combinedObj.integrated[[i]] <- NULL
}

combinedObj.integrated = ScaleData(combinedObj.integrated, verbose = FALSE)
combinedObj.integrated = RunPCA(combinedObj.integrated, verbose = FALSE)
combinedObj.integrated = RunTSNE(combinedObj.integrated, reduction.name = 'rna.tsne', reduction.key = 'rnaTSNE_', dims = 1:30)
combinedObj.integrated <- RunUMAP(combinedObj.integrated, reduction = 'pca', dims = 1:30, assay = 'RNA', 
                                reduction.name = 'rna.umap', reduction.key = 'rnaUMAP_')

combinedObj.integrated <- FindNeighbors(combinedObj.integrated, reduction = "pca", dims = 1:30)
plan("sequential")
combinedObj.integrated <- FindClusters(combinedObj.integrated)



## ----RNA UMAP Plot, echo=FALSE, warning=FALSE, message=FALSE, results="hide"---------
png("UMAP_RNA.png", width=1800, height=1600, res = 300)
DimPlot(combinedObj.integrated, reduction='rna.umap', group.by='integratedRNA_snn_res.0.8', label = TRUE) + ggtitle("RNA")
dev.off()

cols = ceiling(sqrt(length(unique(combinedObj.integrated$orig.ident))))
rows = length(unique(combinedObj.integrated$orig.ident)) / cols
png("UMAP_RNA_Split_Raster.png", width=(1300*cols)+300, height=(1200*rows)+300, res = 300)
print(DimPlot(combinedObj.integrated, reduction='rna.umap', group.by='integratedRNA_snn_res.0.8', split.by='orig.ident', ncol= cols, label = TRUE, raster=TRUE) + ggtitle("RNA"))
dev.off()


if ("BEST" %in% colnames(combinedObj.integrated@meta.data)) {
  png("UMAP_RNA_Demux.png", width=1800, height=1600, res = 300)
  print(DimPlot(combinedObj.integrated, reduction='rna.umap', group.by='BEST', label = TRUE) + ggtitle("RNA"))
  dev.off()
  
  cols = ceiling(sqrt(length(unique(combinedObj.integrated$BEST))))
  rows = length(unique(combinedObj.integrated$BEST)) / cols
  png("UMAP_RNA_Demux_Split_Raster.png", width=(1300*cols)+300, height=(1200*rows)+300, res = 300)
  print(DimPlot(combinedObj.integrated, reduction='rna.umap', split.by='BEST', ncol= cols, label = TRUE, raster=TRUE) + ggtitle("RNA"))
  dev.off()
}

## ----Normalize HTO Data, warning=FALSE, message=FALSE, eval=hashtag------------------
if ('hash.ID' %in% colnames(combinedObj.integrated@meta.data)) {

  ## ----HTO Ridge Plots, results='show', eval=hashtag, fig.width=12, fig.height=9-------
  png('HTO_Ridge_Plot.png', units='in', width=12, height=9, res=300)
  for (i in seq(1,length(rownames(combinedObj.integrated[["HTO"]])), by=25)) {
    print(RidgePlot(combinedObj.integrated, sort(rownames(combinedObj.integrated[['HTO']]))[i:min(i+24,length(rownames(combinedObj.integrated[['HTO']])))], assay="HTO", ncol=min(5, ceiling(sqrt(length(rownames(combinedObj.integrated[['HTO']]))-(i-1)))), group.by='hash.ID'))
  }
  dev.off()
 
  png("UMAP_RNA_HTO.png", width=1800, height=1600, res = 300)
  DimPlot(combinedObj.integrated, reduction='rna.umap', group.by='hash.ID', label = TRUE) + ggtitle("RNA")
  dev.off()
}


################################################################
### B. Performing sample integration on combined object CITESeq assay 
DefaultAssay(object = combinedObj.integrated) <- "ADT"

combinedObj.integrated = NormalizeData(combinedObj.integrated, assay= "ADT", normalization.method = 'CLR', margin = 2, verbose = FALSE)
reference.features = rownames(combinedObj.integrated[["ADT"]])
combinedObj.integrated = ScaleData(combinedObj.integrated, assay= "ADT", verbose = FALSE)
combinedObj.integrated = RunPCA(combinedObj.integrated, assay= "ADT", reduction.name = 'apca', features = reference.features, verbose = FALSE)
combinedObj.integrated <- FindNeighbors(combinedObj.integrated, dims = 1:min(length(rownames(combinedObj.integrated[['ADT']]))-1, 20), reduction = "apca")
combinedObj.integrated <- FindClusters(combinedObj.integrated, graph.name = "ADT_snn", algorithm = 3, verbose = FALSE)
combinedObj.integrated <- RunUMAP(combinedObj.integrated, reduction = 'apca', dims = 1:30, assay = 'ADT', 
                                reduction.name = 'adt.umap', reduction.key = 'adtUMAP_')

png("UMAP_ADT.png", width=1800, height=1600, res = 300)
DimPlot(combinedObj.integrated, reduction='adt.umap', group.by='ADT_snn_res.0.8', label = TRUE) + ggtitle("ADT")
dev.off()

cols = ceiling(sqrt(length(unique(combinedObj.integrated$orig.ident))))
rows = length(unique(combinedObj.integrated$orig.ident)) / cols
png("UMAP_ADT_Split_Raster.png", width=(1300*cols)+300, height=(1200*rows)+300, res = 300)
print(DimPlot(combinedObj.integrated, reduction='adt.umap', split.by='orig.ident', ncol= cols, label = TRUE, raster=TRUE) + ggtitle("ADT"))
dev.off()


################################################################
### C. Multimodal integration for RNA and proteins

future.seed = NULL

DefaultAssay(object = combinedObj.integrated) <- "integratedRNA"

multimode.integrated <- combinedObj.integrated

# Identify multimodal neighbors. These will be stored in the neighbors slot, and can be accessed using bm[['weighted.nn']]
# The WNN graph can be accessed at bm[["wknn"]], and the SNN graph used for clustering at bm[["wsnn"]]
# Cell-specific modality weights can be accessed at bm$RNA.weight

### Currently commenting out because it is slow and not always useful

# multimode.integrated <- FindMultiModalNeighbors(
#   combinedObj.integrated, reduction.list = list("pca", "apca"), 
#   dims.list = list(1:30, 1:30), modality.weight.name = c("intRNA.weight", "intADT.weight"))
# 
# multimode.integrated <- RunUMAP(multimode.integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
# 
# multimode.integrated <- FindClusters(multimode.integrated, graph.name = "wsnn", verbose = FALSE)
# 
# 
# 
# p1 <- DimPlot(multimode.integrated, reduction = 'wnn.umap', label = TRUE, repel = TRUE, label.size = 2.5) #+ NoLegend()
# 
# png("UMAP_MultiModal_WNN.png", width=1800, height=1600, res = 300)
# p1 + ggtitle("Multi-Modal Weighted Graph")
# dev.off()
# 
# cols = ceiling(sqrt(length(unique(multimode.integrated$orig.ident))))
# rows = length(unique(multimode.integrated$orig.ident)) / cols
# png("UMAP_MultiModal_WNN_Split_Raster.png", width=(1300*cols)+300, height=(1200*rows)+300, res = 300)
# print(DimPlot(multimode.integrated, reduction='wnn.umap', split.by='orig.ident', ncol= cols, label = TRUE, raster=TRUE) + ggtitle("Multi-Modal WNN"))
# dev.off()
# 
# 
# ## ----Multi-Modal Modality Weights, echo=FALSE, warning=FALSE, message=FALSE, results='hide', fig.width=10, fig.height=7----
# png("MultiModal_WNN_Weights.png", width=3600, height=2400, res = 300)
# (VlnPlot(multimode.integrated, features = "SCT.weight", group.by = 'wsnn_res.0.8', pt.size = 0.1) + NoLegend()) / (VlnPlot(multimode.integrated, features = "ADT.weight", group.by = 'wsnn_res.0.8', pt.size = 0.1) + NoLegend())
# dev.off()
# 
# 
# p3 <- DimPlot(multimode.integrated, reduction = 'rna.umap', label = TRUE, 
#               repel = TRUE, label.size = 2.5) + NoLegend() + ggtitle("RNA")
# p4 <- DimPlot(multimode.integrated, reduction = 'adt.umap', label = TRUE, 
#               repel = TRUE, label.size = 2.5)  + ggtitle("ADT") #+ NoLegend()
# png("UMAP_MultiModal_RNA_ADT.png", width=3600, height=1600, res = 300)
# p3 + p4
# dev.off()

saveRDS(multimode.integrated, "multimode.integrated.rds")


## ----Prepare for Heatmap, include=FALSE----------------------------------------------
require(scales)
identities <- levels(multimode.integrated$integratedRNA_snn_res.0.8)
DefaultAssay(multimode.integrated) <- 'integratedRNA'
Idents(multimode.integrated) <- multimode.integrated$integratedRNA_snn_res.0.8

suppressWarnings(dir.create('DEG'))


## ----DEG Calculation, message=FALSE, results='asis'----------------------------------

sample <- 'SeuratAggregate'
# 1 cluster vs all other clusters
for (cluster in c(0:(length(identities) - 1))) {
  try({
    Markers <- FindMarkers(multimode.integrated,ident.1=cluster,test.use="MAST",only.pos = T)
    if (length(row.names(Markers) != 0)){
      Markers$sub <- Markers$pct.1-Markers$pct.2
      Markershead <- head(Markers, n=20)
      
      filename <- paste0('DEG/', sample, "_", cluster, "_markers.txt")
      write.table(cbind(Genes = rownames(Markers), Markers), filename, sep="\t", quote=FALSE, row.names=FALSE)
      
      Markers12 <- head(rownames(Markershead), n=12)
      filename <- paste0('DEG/', sample, "_", cluster,"_Features.png", sep="")
      png(filename , width=3200, height=3200, res = 300)
      print(FeaturePlot(multimode.integrated, features = Markers12, ncol=3))
      dev.off()
    }
  })
}
print('DONE')
