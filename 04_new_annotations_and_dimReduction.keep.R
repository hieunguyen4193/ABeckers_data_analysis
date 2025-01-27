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

gene.list <- list(
  cDC2B = c("B2M", "FTH1", "LYZ", "S100A9", "S100A8", "VCAN",  "ANXA1", "CD14",  "NOTCH2", "CLEC7A", "ITGAM", "KLF4",  "SIRPA", "CLEC4A", 
            "MT2A", "CD52", "OGFRL1", "CEBPD", "LGALS2", "CEBPD",  "CRYBG1", "JUNB", "MAP3K8", "ATP1B3", "GPR183", "TNFAIP3", "SLC38A1", "RGS2", 
            "RTN1", "IL1B", "TEX14", "CCNH",  "GADD45B",  "CD1C",  "CLEC10A",   "FCER1A", "ENHO",  "CD2", "CD3E", "LTB"),
  cDC2A = c("ITGA4", "VAV3", "ZCCHC7", "AFF1", "SSBP2", "CDK6", "AUTS2", "RUNX2", "BCL11A", "TCF4", "IGKC",  "ITM2C",   "LTB", "P2RY14", "NUCB2"),
  DC3 = c("S100A12", "LYST", "MBNL1", "LRMDA",  "PLCB1", "ARHGAP26", "RBM47", "ZSWIM6", "GNAQ",  "DIAPH2", "SIPA1L1", "ASAP1",  "FNDC3B",  "SYK",   "GPCPD1","AGTPBP1", "ATXN1",     "FRY",  "SLC39A11", "SGMS1", "SNTB1"),
  DC1 = c("CD74", "S100A10", "CST3", "HLA-DRA",  "GSTP1", "CRIP1",  "HLA-DPA1", "HLA-DPB1", "SNX3", "UCP2", "EEF1B2", "PSMB9", "CD63", "ID2", "UPF2", "RGS10", "FNBP1", "CPVL", "PLEK", "CD83", "HLA-DQB1", "HLA-DRB5", "HLA-DQA1", "HLA-DMA","RAB32", "DAPP1", "SLAMF7", "TFRC", "SIPA1L3", "VOPP1", "CAMK2D",  "SHTN1", "KAT2B",  "SECISBP2L",  "WDFY4",  "CPNE3", "ITGAE", "GPR137B", "BATF3", "THBD", "CLEC9A", "CD226", "CD59", "C1orf54", "CCSER1", "DUSP2", "PLEKHA5", "DDIT4", "SLAMF8", "NDRG2", "IDO1", "CLNK", "ZNF366", "XCR1", "DPP4", "BTLA", "DPP4", "CADM1",   "DNASE1L3"),
  pDC = c("TCF4", "LILRA4", "IL3RA", "CLEC4C", "TLR9", "NRP1", "FCER1A", "CD2", "CD81", "BCL11A"),
  preDC = c("PTPRC", "ZEB2", "SPI1", "IKZF1", "KLF4",  "SEMA4D", "RELB", "IRF5", "CSF1R", "CEBPA", "IRF8",  "IRF7", "IL3RA", "FLT3", "ZBTB46", "CLEC4C", "IRF4", "SPINK2", "IL7R", "BCL2")
)


test.gene.list <- c(
  "PTPRC", "ZEB2", "SPI1", "IKZF1", "KLF4",  "SEMA4D", "RELB", "IRF5", "CSF1R", "CEBPA", "IRF8",  "IRF7", "IL3RA", "FLT3", "ZBTB46", "CLEC4C", "IRF4", "SPINK2", "IL7R", "BCL2"
)
pca_reduction_name <- "PCA_subset"
tsne_reduction_name <- "TSNE_subset"
my_random_seed <- 42
s.obj <- subset(s.obj, seurat_clusters == 2)

DefaultAssay(s.obj) <- "RNA"

# s.obj <- ScaleData(s.obj)
# s.obj <- NormalizeData(s.obj)
# s.obj <- FindVariableFeatures(s.obj)

s.obj <- RunPCA(s.obj, 
                npcs = num.PCA, 
                verbose = FALSE, 
                reduction.name=pca_reduction_name,
                features = test.gene.list)

get.real.PCA.dim <- ncol(s.obj@reductions$PCA_subset@cell.embeddings)
s.obj <- RunTSNE(s.obj, reduction = pca_reduction_name, 
                 dims = 1:get.real.PCA.dim, reduction.name=tsne_reduction_name,
                 seed.use = my_random_seed, check_duplicates = FALSE)
# clustering 
s.obj <- FindNeighbors(s.obj, reduction = pca_reduction_name, dims = 1:get.real.PCA.dim)
s.obj <- FindClusters(s.obj, resolution = cluster.resolution, random.seed = 0)

DimPlot(object = s.obj, reduction = "TSNE_subset", label = TRUE, label.box = TRUE)

FeaturePlot(object = s.obj, 
            features = head(test.gene.list, 9), 
            ncol = 3, 
            reduction = "TSNE_subset",
            label = TRUE)
