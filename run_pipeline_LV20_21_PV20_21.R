gc()
rm(list = ls())
my_random_seed <- 42

set.seed(my_random_seed)
# __________VDJ DATA ANYLYSIS PIPELINE__________
PROJECT <- "ABeckers_data"
batch.id <- "batch1"

# __________GEX DATA ANALYSIS PIPELINE__________
path.to.pipeline.src <- "/home/hieunguyen/CRC1382/src_2023/src_pipeline/scRNA_GEX_pipeline"
path2src <- file.path(path.to.pipeline.src, "processes_src")

source(file.path(path2src, "import_libraries.R"))
source(file.path(path.to.pipeline.src, "scRNA_GEX_pipeline.R"))
path.to.storage <- "/media/hieunguyen/HD01/storage"
path2input <- file.path(path.to.storage, PROJECT, batch.id)
outdir <- "/home/hieunguyen/CRC1382/outdir"

# output.version <- "default"
# output.version <- "v0.1"
cluster.resolution <- 1

for (output.version in c("default", "v0.1")){
  source("/home/hieunguyen/CRC1382/src_2023/ABeckers/config.R")
  
  path.to.output <- file.path(outdir, PROJECT, output.version, batch.id, sprintf("res_%s", cluster.resolution))
  dir.create(path.to.output, showWarnings = FALSE, recursive = TRUE)
  
  # _____stage lst for single sample_____
  stage_lst <- list()
  
  stage_lst = c(
    LV20      =   "LV",
    LV21      =   "LV",
    PV20      =   "PV",
    PV21      =   "PV"
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
  
  num.PCA <- 30
  num.PC.used.in.UMAP <- 30
  num.PC.used.in.Clustering <- 30
  
  num.dim.integration <- 30
  num.dim.cluster <- 30
  
  filtered.barcodes <- NULL
  
  s.obj <- run_pipeline_GEX(path2src=path2src,
                            path2input=path2input,
                            path.to.logfile.dir=file.path(path.to.output, "logs"),
                            stage_lst=stage_lst,
                            path.to.10X.doublet.estimation=path.to.10X.doublet.estimation,
                            MINCELLS=MINCELLS,
                            MINGENES=MINGENES,
                            PROJECT=PROJECT,
                            remove_doublet=remove_doublet,
                            save.RDS=save.RDS,
                            path.to.output=path.to.output,
                            rerun=rerun,
                            DE.test="wilcox",
                            num.PCA=num.PCA,
                            num.PC.used.in.UMAP=num.PC.used.in.UMAP,
                            num.PC.used.in.Clustering=num.PC.used.in.Clustering,
                            use.sctransform=FALSE,
                            filtered.barcodes=filtered.barcodes,
                            filter.thresholds=filter.thresholds,
                            path.to.anno.contigs=path.to.anno.contigs,
                            path.to.count.clonaltype=path.to.count.clonaltype,
                            input.method = "CITESEQ",
                            my_random_seed = my_random_seed,
                            cluster.resolution = cluster.resolution, 
                            num.dim.integration = num.dim.integration,
                            sw = sw, with.TSNE = TRUE)
  
  #### ALWAYS REMEMBER TO SAVE SESSIONINFO !!!!!!
  writeLines(capture.output(sessionInfo()), file.path(path.to.output, sprintf("%s_sessionInfo.txt", PROJECT)))
}


