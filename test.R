library(Seurat)

tmp.old <- readRDS("/media/hieunguyen/HD01/outdir/CRC1382/AGBerres/220311_Berres_MiedIII_scCITEseq/s8_output/AGBerres_4_samples_LV_PV.output.s8.rds")
tmp.new <- readRDS("/home/hieunguyen/CRC1382/outdir/ABeckers_data/batch1/s8_output/ABeckers_data.output.s8.rds")

# count.mat.old <- GetAssayData(object = tmp.old, assay = "integrated",slot = "data")
count.mat.new <- GetAssayData(object = tmp.new, assay = "integrated",slot = "data")

num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
cluster.resolution <- 0.8
num.dim.integration <- 30
num.dim.cluster <- 30
cluster.resolution <- 0.5
inte_pca_reduction_name = "INTE_PCA"
inte_umap_reduction_name = "INTE_UMAP"

s.obj <- tmp.new
DefaultAssay(s.obj) <- "RNA"

data.list <- SplitObject(s.obj, split.by = "name")

data.list <- lapply(X = data.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nåfeatures = 2000)})    

anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:num.dim.integration, scale=F)## THIS IS CCA DIMENSIONS

s.obj_inte <- IntegrateData(anchorset = anchors, dims = 1:num.dim.integration) ## THIS IS PCA DIMENSION

## keep the order of integration obj
s.obj_inte <- s.obj_inte[, colnames(s.obj)]

s.obj[['integrated']] <- s.obj_inte[['integrated']]

s.obj@commands <- c(s.obj@commands, s.obj_inte@commands)

s.obj@tools <- c(s.obj@tools, s.obj_inte@tools)

DefaultAssay(s.obj) <- "integrated"

s.obj <- ScaleData(s.obj, verbose = FALSE, features = row.names(s.obj))

s.obj <- FindVariableFeatures(s.obj)
pca.genes <- setdiff(VariableFeatures(s.obj), genes.to.not.run.PCA)
print("####################################################################")
print(sprintf("Running PCA with %s genes, after removing %s genes", length(pca.genes), length(genes.to.not.run.PCA)))
print("####################################################################")
s.obj <- RunPCA(s.obj, npcs = num.PCA, verbose = FALSE, reduction.name=inte_pca_reduction_name, 
                features = pca.genes)

s.obj <- RunUMAP(s.obj, reduction = inte_pca_reduction_name, dims = 1:num.PC.used.in.UMAP, reduction.name=inte_umap_reduction_name, 
                 seed.use = my_random_seed, umap.method = umap.method)