gc()
rm(list = ls())

#####----------------------------------------------------------------------#####
##### Install new packages
#####----------------------------------------------------------------------#####
list.of.packages <- c("Seurat",
                      "SingleCellExperiment",
                      "optparse", 
                      "comprehenr", 
                      "tidyverse", 
                      "ggplot2", 
                      "SoupX",
                      "comprehenr",
                      "DoubletFinder",
                      "vroom",
                      "hash",
                      "DT",
                      "janitor",
                      "knitr",
                      "circlize",
                      "formattable",
                      "htmlwidgets",
                      "plotly",
                      "stringr"
)

bioc.packages <- c("celda", 
                   "BiocSingular", 
                   "PCAtools", 
                   "SingleCellExperiment",
                   "scRepertoire")

# Check if packages are installed 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
new.bioc.packages <- bioc.packages[!(bioc.packages %in% installed.packages()[,"Package"])]

# Install new packages 
if(length(new.packages)) install.packages(new.packages)
if(length(new.bioc.packages)) BiocManager::install(new.bioc.packages)

# Import all packages 
package_loading_Status <- lapply(list.of.packages, 
                                 require, 
                                 character.only = TRUE)

package_loading_Status_bioc <- lapply(bioc.packages, 
                                      require, 
                                      character.only = TRUE)

#####----------------------------------------------------------------------#####
##### MAIN FUNCTIONS
#####----------------------------------------------------------------------#####

# path to the RDS file
path.to.s.obj <- "/home/hieunguyen/CRC1382/outdir/merge_all_6_samples_test_function.output.s8.rds"

# path to save output svg file
path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/tmp"
dir.create(path.to.save.output, showWarnings = FALSE, recursive = TRUE)

s.obj <- readRDS(path.to.s.obj)

Idents(s.obj) <- "celltype" # <<<< CHANGE HERE to use cluster number.
# Idents(s.obj) <- "seurat_clusters" # <<<<< CHANGE HERE to switch to cell annotation

# set default assay =  RNA
DefaultAssay(s.obj) <- "RNA"
reduction.name <- "INTE_TSNE"
# genes to show on the UMAP, violin plot
plot.genes <- c("ID2", "BATF3", "ZEB2")

##### UMAP plot, all samples
umap.plot <- DimPlot(object = s.obj, reduction = reduction.name, label = TRUE, label.box = TRUE, repel = TRUE, pt.size = 1)

##### UMAP plot, splitted by names
umap.plot.names <- DimPlot(object = s.obj, 
                           reduction = reduction.name,
                           label = TRUE,
                           label.box = FALSE, 
                           repel = TRUE, 
                           split.by = "name", 
                           ncol = 2,
                           pt.size = 1,
                           label.size = 12)

##### feature plot
feature.plot.sample <- FeaturePlot(object = s.obj, reduction = reduction.name, label = TRUE, repel = TRUE, ncol = 2, pt.size = 0, features = plot.genes, slot = "data")

##### violin plot
violin.plot <- VlnPlot(object = s.obj, features = plot.genes, ncol = 2, pt.size = 0, slot = "data")  # set pt.size = 0 to remove points on violin plots

##### save function
choose_your_filename <- "feature_plot_01.svg"
ggsave(plot = feature.plot.sample, filename = choose_your_filename, path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)

choose_your_filename <- "violin_plot_01.svg"
ggsave(plot = violin.plot, filename = choose_your_filename, path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)

choose_your_filename <- "heatmap_plot_01.svg"
hm <- DoHeatmap(object = s.obj, features = plot.genes, assay = "RNA", draw.lines = TRUE, slot = "data") +
  scale_fill_gradient(low = "white", high = "red")
ggsave(plot = hm, filename = choose_your_filename, path = path.to.save.output, device = "svg", width = 20, height = 15, dpi = 300)
