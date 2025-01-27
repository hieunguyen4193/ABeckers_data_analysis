##### clean up #####
gc()
rm(list = ls())

path.to.project.src <- "/home/hieunguyen/CRC1382/src_2023/ABeckers"
source(file.path(path.to.project.src, "00_import_libraries.R"))
source(file.path(path.to.project.src, "00_helper_functions.R"))

path.to.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline/processes_src"
source(file.path(path.to.src, "s8_integration_and_clustering.R"))

#####----------------------------------------------------------------------#####
# CONFIGURATIONS AND PREPRATIONS
#####----------------------------------------------------------------------#####
PROJECT <- "ABeckers_data"
outdir <- "/home/hieunguyen/CRC1382/outdir"

output.version <- "default"
cluster.resolution <- 0.8
path.to.main.input <- file.path(outdir, PROJECT)
path.to.main.output <- file.path(path.to.main.input, output.version, "integrated_data", sprintf("res_%s", cluster.resolution))

path.to.03.output <- file.path(path.to.main.output, "03_output")
path.to.04.output <- file.path(path.to.main.output, "04_output")
dir.create(path.to.04.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(file.path(path.to.03.output, 
                           "s8_output", 
                           "merge_all_6_samples_test_function.output.s8.rds"))

##### delete cluster B cell
s.obj <- subset(s.obj, celltype != "B-cell")

num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
num.dim.integration <- 30
num.dim.cluster <- 30
chosen.seed <- 42

gene.list <- c( "CRIP1", "CD83", "HLA-DQB1", "HLA-DRB5", "HLA-DQA1", "JUNB", "VAV3", "RTM1", "TEX14", "SSBP2", "LRMDA", "PTCB1", "FCAR")

DefaultAssay(s.obj) <- "RNA"
s.obj <- AddModuleScore(object = s.obj, features = list(DiffGene = intersect(gene.list, row.names(s.obj))), name = sprintf("%s_", "DiffGene"), ctrl = 50)

feature.plot.module.scores <- FeaturePlot(object = s.obj, 
                                          reduction = "INTE_TSNE", 
                                          label = TRUE, 
                                          features = "DiffGene_1", 
                                          pt.size = 0.5) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")

VlnPlot(object = s.obj, features = c("DiffGene_1"))

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") 
meta.data$seurat_clusters <- as.numeric(meta.data$seurat_clusters)

meta.data$seurat_clusters_new <- unlist(mapply(function(x, y){
  if (x %in% c(4, 6, 9)){
    subset.metadata <- subset(meta.data, meta.data$seurat_clusters %in% c(x))
    thres <- quantile(subset.metadata$DiffGene_1)
    median.diff.gene <- thres[["50%"]]
    if (y >= median.diff.gene){
      output <- sprintf("%sA", x)
    } else {
      output <- sprintf("%sB", x)
    }
  } else if (x %in% c(4, 6, 9) == FALSE) {
    output <- sprintf("%s", x)
  }
  return(output)
}, meta.data$seurat_clusters, meta.data$DiffGene_1 ))

cluster.convert <- list(
  `1` = "cMo/cMo_M1",
  `2` = "cDC1",
  `3` = "cDC3",
  `4A` = "cMo/cMo_M1", 
  `4B` = "iMo_M2",
  `5` = "pDC precursor",
  `6A`= "cMo/cMo_M1",
  `6B` = "iMo_M2",
  `7` = "cMo/cMo_M1",
  `8` = "pDC",
  `9A`= "cMo/cMo_M1", 
  `9B` = "iMo_M2",
  `10` = "iMo_M2",
  `11` = "cDC precursor",
  `12` = "cDC2",
  `13`= "cMo/cMo_M1"
)

meta.data <- meta.data %>% rowwise() %>% 
  mutate(cluster.name = cluster.convert[[seurat_clusters_new]])

meta.data <- meta.data %>% column_to_rownames("barcode")
meta.data <- meta.data[row.names(s.obj@meta.data), ]

s.obj <- AddMetaData(object = s.obj, metadata = meta.data$cluster.name, col.name = "new_annotation")
new.tsne.plot <- DimPlot(object = s.obj, reduction = "INTE_TSNE", label = TRUE, group.by = "new_annotation", label.box =  TRUE)

saveRDS(s.obj, file.path(path.to.04.output, "merge_all_6_samples_test_function.output.s8.newAnnotation.rds"))

# pca_reduction_name <- "PCA_subset"
# tsne_reduction_name <- "TSNE_subset"
# my_random_seed <- 42
# s.obj <- RunPCA(s.obj, 
#                 npcs = num.PCA, 
#                 verbose = FALSE, 
#                 reduction.name=pca_reduction_name,
#                 features = gene.list)
# 
# get.real.PCA.dim <- ncol(s.obj@reductions$PCA_subset@cell.embeddings)
# s.obj <- RunTSNE(s.obj, reduction = pca_reduction_name, 
#                  dims = 1:get.real.PCA.dim, reduction.name=tsne_reduction_name,
#                  seed.use = my_random_seed, check_duplicates = FALSE)
# 
# # clustering 
# s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:get.real.PCA.dim)
# s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)
# 
# DimPlot(object = s.obj, reduction = "TSNE_subset", label = TRUE, label.box = TRUE)
# 
# FeaturePlot(object = s.obj, 
#             features = head(gene.list, 9), 
#             ncol = 3, 
#             reduction = "TSNE_subset",
#             label = TRUE)