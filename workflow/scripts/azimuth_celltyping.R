library(Seurat)
library(Azimuth)
library(SeuratData)
library(ggplot2)


args <- commandArgs(TRUE);
workdir <- args[1]
rds_path <- args[2]
genome <- args[3]
ref <- args[4]
azimuth_ref <- args[5]


seur <- readRDS(rds_path)
if (ref == 'pbmc' & genome == 'hg38' & (!is.na(azimuth_ref) & azimuth_ref != '')){
  # reference = "/data/OpenOmics/references/chicyte/hg38/Azimuth/pbmcref/"
  seur <- RunAzimuth(seur, reference = azimuth_ref)
}

groups <- grep('score', grep('predicted.celltype', names(seur@meta.data), value=TRUE), value=TRUE, invert=TRUE)

## ----Create Working Directory, include=FALSE-----------------------------------------
dir.create(workdir, recursive=TRUE)


## ----Change Working Directory, include=FALSE-----------------------------------------
setwd(workdir)


#p2 <- DimPlot(pbmcsca, group.by = "Method")
#p1 + p2

for (group in groups) {
  level <- strsplit(group, '.', fixed=TRUE)[[1]][3]
  width <- 1200+400*round(dim(unique(seur[[group]]))[[1]]/12)
  png(paste0('UMAP_RNA_CellType_', level, '.png'), width=width, height=1600, res = 300)
  print(DimPlot(seur, group.by = group, label = TRUE, label.size = 3, repel = TRUE, raster = TRUE) + ggtitle("PBMC"))
  dev.off()

  png(paste0('CellScoreViolinPlot_', level, '.png'), width=1200+800*round(dim(unique(seur[[group]]))[[1]]/12), height=1600, res = 300)
  print(VlnPlot(seur, features=paste0(group, '.score'), group.by=group))
  dev.off()

  write.table(table(seur[[group]]), paste0("CellTypeCount_", level, ".csv"), quote=FALSE, row.names=FALSE, col.names=FALSE, sep=',')

  if (dim(unique(seur[['orig.ident']]))[[1]] > 1) {
    png(paste0('UMAP_RNA_CellType_', level, '_Split.png'), width=width+400*dim(unique(seur[['orig.ident']]))[[1]], height=1600, res = 300)
    print(DimPlot(seur, group.by = group, split.by = 'orig.ident', label = TRUE, label.size = 3, repel = TRUE, raster = TRUE) + ggtitle("PBMC"))
    dev.off()
  }
}



saveRDS(seur, 'azimuth_prediction.rds')
