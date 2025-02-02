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

Idents(s.obj) <- "celltype"
feature.plot.module.scores <- FeaturePlot(object = s.obj, 
                                          reduction = "INTE_TSNE", 
                                          label = TRUE, 
                                          features = "DiffGene_1", 
                                          pt.size = 0.5) &
  scale_color_gradient(low = "lightgray", high = "#FF0000", na.value = "lightgray")

to.recluster <- c(
  "MoMac_0", "MoMac_2", "MoMac_7"
)

meta.data <- s.obj@meta.data %>% rownames_to_column("barcode") 

convert.celltype <- list(
  `Momac_1` = "cMo/M2",
  `DC3` = "cDC3",
  `DC2` = "cDC2",
  `NCM` = "cMo/M2",
  `DC1`	= "cDC1",
  `MoMac_6` = "cMo/M2",
  `pDC` = "pDC",
  `MDP 16+` = "pDC precursor",
  `MDP` = "precursor",
  `MoMac_12` = "cMo/M1",
  `CDP` = "cDC precursor"
  )

tex14.exprs <- GetAssayData(object = s.obj, slot = "counts", assay = "RNA")["TEX14", ]
tex14.pos <- tex14.exprs[tex14.exprs != 0] %>% names()
tex14.zero <- tex14.exprs[tex14.exprs == 0] %>% names()

new.annotations <- unlist(mapply(function(x, y, z){
  if (x %in% to.recluster){
    subset.metadata <- subset(meta.data, meta.data$celltype == x)
    thres <- quantile(subset.metadata$DiffGene_1)
    median.diff.gene <- thres[["50%"]]
    if (y >= median.diff.gene){
      output <- "cMo/M2"
    } else {
      output <- "cMo/M1"
    }
  } else if (x %in% to.recluster == FALSE && x != "DC2") {
    output <- convert.celltype[[x]]
  } else if (x == "DC2"){
    if (z %in% tex14.pos){
      output <- "cDC2A"
    } else {
      output <- "cDC2B"
    }
  }
  return(output)
}, meta.data$celltype, meta.data$DiffGene_1, meta.data$barcode))

meta.data$new_annotation <- new.annotations
meta.data <- meta.data %>% column_to_rownames("barcode")

meta.data <- meta.data[rownames(s.obj@meta.data), ]
s.obj <- AddMetaData(object = s.obj, metadata = meta.data$new_annotation, col.name = "new_annotation")
# DimPlot(object = s.obj, reduction = "INTE_TSNE", label = TRUE, label.box = TRUE, repel = TRUE, group.by = "new_annotation")

saveRDS(s.obj, file.path(path.to.04.output, "merge_all_6_samples_test_function.output.s8.newAnnotation_20250129.rds"))
