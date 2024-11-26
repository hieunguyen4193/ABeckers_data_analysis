

gc()
rm(list = ls())
my_random_seed <- 42

set.seed(my_random_seed)
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "ABeckers_data"
batch.id <- "batch2"

num.PCA <- 30
num.PC.used.in.UMAP <- 30
num.PC.used.in.Clustering <- 30
# cluster.resolution <- 0.8
cluster.resolution <- 1
num.dim.integration <- 30
num.dim.cluster <- 30
input.method <- "CITESEQ"

path.to.storage <- "/media/hieunguyen/HD01/storage"
outdir <- "/home/hieunguyen/CRC1382/outdir"
source("/home/hieunguyen/CRC1382/src_2023/ABeckers/config.R")

# output.version <- "v0.1"
for (output.version in c("v0.1", "default")){
  # __________GEX DATA ANALYSIS PIPELINE__________
  path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
  path2src <- file.path(path.to.pipeline.src, "processes_src")
  
  source(file.path(path2src, "import_libraries.R"))
  source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))
  
  path2input <- file.path(path.to.storage, PROJECT, batch.id)
  path.to.output <- file.path(outdir, PROJECT, output.version, batch.id, sprintf("res_%s", cluster.resolution))
  
  # _____stage lst for single sample_____
  stage_lst <- list()
  
  stage_lst = c(
    LV73      =   "LV",
    PV73      =   "PV"
  )
  
  MINCELLS  <- 5
  MINGENES  <- 50
  
  save.RDS <- list(s1 = TRUE,
                   s2 = TRUE,
                   s3 = TRUE,
                   s4 = TRUE,
                   s5 = TRUE,
                   s6 = TRUE,
                   s7 = FALSE,
                   s8 = TRUE,
                   s8a = TRUE,
                   s9 = TRUE)
  
  sw <- list(s1 = "on",
             s2 = "on",
             s3 = "on",
             s4 = "on",
             s5 = "on",
             s6 = "on",
             s7 = "off",
             s8 = "on",
             s8a = "off",
             s9 = "on")
  
  rerun <- list(s1 = FALSE, 
                s2 = FALSE,
                s3 = FALSE,
                s4 = FALSE,
                s5 = FALSE,
                s6 = FALSE,
                s7 = FALSE,
                s8 = FALSE,
                s8a = FALSE,
                s9 = FALSE)
  
  
  filter.thresholds <- all.config.params[[output.version]]
  
  print("Using the following configurations:")
  for (i in names(filter.thresholds)){
    print(sprintf("Param %s: %s", i, filter.thresholds[[i]]))
  }
  
  remove_doublet <- FALSE
  path.to.10X.doublet.estimation <- "/media/hieunguyen/HD01/storage/DoubletEstimation10X.csv"
  
  filtered.barcodes <- NULL
  dir.create(file.path(path.to.output, sprintf("res_%s", cluster.resolution)), showWarnings = FALSE, recursive = TRUE)
  s.obj <- run_pipeline_GEX(path2src = path2src,
                            path2input = path2input,
                            path.to.logfile.dir = file.path(path.to.output, "logs"),
                            stage_lst = stage_lst,
                            path.to.10X.doublet.estimation = path.to.10X.doublet.estimation,
                            MINCELLS = MINCELLS,
                            MINGENES = MINGENES,
                            PROJECT = PROJECT,
                            remove_doublet = remove_doublet,
                            save.RDS = save.RDS,
                            path.to.output = path.to.output,
                            rerun = rerun, 
                            DE.test = "wilcox",
                            num.PCA = num.PCA,
                            num.PC.used.in.UMAP = num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering = num.PC.used.in.Clustering,
                            use.sctransform = FALSE,
                            filtered.barcodes = filtered.barcodes,
                            filter.thresholds = filter.thresholds,
                            input.method = input.method,
                            cluster.resolution = cluster.resolution,
                            num.dim.integration = num.dim.integration,
                            inte_pca_reduction_name = "INTE_PCA",
                            inte_umap_reduction_name = "INTE_UMAP",
                            with.VDJ = FALSE, 
                            sw = sw)
  
  #### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
  writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
  
}

