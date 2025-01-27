gc()
rm(list = ls())

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)

path.to.s.obj <- "/home/hieunguyen/CRC1382/outdir/ABeckers_data/default/integrated_data/res_0.8/04_output/merge_all_6_samples_test_function.output.s8.newAnnotation.rds"
path.to.save.output <- "/home/hieunguyen/CRC1382/outdir/ABeckers_data/default/integrated_data/res_0.8/04_output/tmp_output"

dir.create(path.to.save.output)

s.obj <- readRDS(path.to.s.obj)
DefaultAssay(s.obj) <- "RNA"

Idents(s.obj) <- "new_annotation"
s.obj.LV <- subset(s.obj, stage == "LV")
s.obj.PV <- subset(s.obj, stage == "PV")

cluster_of_interest <- c("cDC1",
                         "cMo/cMo_M1",
                         "cDC3",
                         "iMo_M2",
                         "pDC precursor",
                         "pDC",
                         "cDC precursor")
gene_of_interest <- list()
gene_of_interest[["Maturation"]] <- c("CD40", "CD86", "RELB", "CD83")
gene_of_interest[["Regulatory"]] <- c("CD274", "CD200", "FAS", "ALDH1A2", "SOCS1", "SOCS2", "BIRC3", "MARCKSL1")
gene_of_interest[["TLRs_and_Adaptors"]] <- c("TICAM1", "TLR9", "TLR8", "TLR7", "TLR6", "TLR5", "TLR4", "TLR2", "TLR1")
gene_of_interest[["Migration"]] <- c("CCR7", "MYO1G", "CXCL16", "ADAM8", "ICAM1", "FSCN1", "MARCKS", "MARCKSL1")

## Internal comparison between clusters in each LV, PV 2-sample 

### LV samples 
DefaultAssay(s.obj) <- "RNA"
for (group in names(gene_of_interest)){
  cat(sprintf("#### %s \n", group))
  lv.df <- data.frame(gene_of_interest[[group]])
  colnames(lv.df) <- c("Gene")
  for (chosen.cluster in cluster_of_interest){
    tmp.s.obj.lv <- subset(s.obj.LV, cells = WhichCells(object = s.obj.LV, ident = chosen.cluster))
    count.matrix <- as.data.frame(GetAssayData(tmp.s.obj.lv, slot = "data"))
    count.matrix <- count.matrix[gene_of_interest[[group]], ]
    countdf <- as.data.frame(rowMeans(count.matrix))
    colnames(countdf) <- c(chosen.cluster)
    countdf <- countdf %>% rownames_to_column("Gene")
    lv.df <- merge(lv.df, countdf, by.x = "Gene", by.y = "Gene")
  }
  lv.df <- lv.df %>% column_to_rownames("Gene")
  lv.df$total <- rowSums(lv.df)
  lv.df <- subset(lv.df, lv.df$total != 0)
  lv.df <- subset(lv.df, select = -c(total))
  input.to.heatmap.lv <- (lv.df - t(rowMeans(lv.df)))/(rowSds(as.matrix(lv.df)))
  
  tmp.heatmap.lv <- input.to.heatmap.lv %>% 
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-c(Gene), names_to = "Cluster", values_to = "Expression") %>%
    mutate(Cluster= fct_relevel(Cluster,colnames(input.to.heatmap.lv))) %>%
    ggplot(aes(x=Cluster, y=Gene, fill=Expression)) + 
    geom_raster() + 
    scale_fill_gradient2(low="#1417f5", high="#f52314", guide="colorbar") +
    theme(text = element_text(size=25), axis.text.x = element_text(angle = 90))
  ggsave(plot = tmp.heatmap.lv, filename = "HEATMAP_LV.svg", path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)
  cat("\n \n")
}

### PV samples 
DefaultAssay(s.obj) <- "RNA"
for (group in names(gene_of_interest)){
  cat(sprintf("#### %s \n", group))
  ##### PREPARE HEATMAP DATA FOR PV SAMPLES
  pv.df <- data.frame(gene_of_interest[[group]])#
  colnames(pv.df) <- c("Gene")
  for (chosen.cluster in cluster_of_interest){
    tmp.s.obj.pv <- subset(s.obj.PV, cells = WhichCells(object = s.obj.PV, ident = chosen.cluster))
    count.matrix <- as.data.frame(GetAssayData(tmp.s.obj.pv, slot = "data"))
    count.matrix <- count.matrix[gene_of_interest[[group]], ]
    countdf <- as.data.frame(rowMeans(count.matrix))
    colnames(countdf) <- c(chosen.cluster)
    countdf <- countdf %>% rownames_to_column("Gene")
    pv.df <- merge(pv.df, countdf, by.x = "Gene", by.y = "Gene")
  }
  pv.df <- pv.df %>% column_to_rownames("Gene")
  pv.df$total <- rowSums(pv.df)
  pv.df <- subset(pv.df, pv.df$total != 0)
  pv.df <- subset(pv.df, select = -c(total))
  input.to.heatmap.pv <- (pv.df - t(rowMeans(pv.df)))/(rowSds(as.matrix(pv.df)))
  
  tmp.heatmap.pv <- input.to.heatmap.pv %>% 
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-c(Gene), names_to = "Cluster", values_to = "Expression") %>%
    mutate(Cluster= fct_relevel(Cluster,colnames(input.to.heatmap.pv))) %>%
    ggplot(aes(x=Cluster, y=Gene, fill=Expression)) + 
    geom_raster() + 
    scale_fill_gradient2(low="#1417f5", high="#f52314", guide="colorbar") +
    theme(text = element_text(size=25), axis.text.x = element_text(angle = 90))
  ggsave(plot = tmp.heatmap.pv, filename = "HEATMAP_PV.svg", path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)
  cat("\n \n")
}

## Cross condition comparisons on heatmaps {.tabset}
DefaultAssay(s.obj) <- "RNA"
for (group in names(gene_of_interest)){
  cat(sprintf("### %s \n", group))
  lv.df <- data.frame(gene_of_interest[[group]])
  colnames(lv.df) <- c("Gene")
  for (chosen.cluster in cluster_of_interest){
    tmp.s.obj.lv <- subset(s.obj.LV, cells = WhichCells(object = s.obj.LV, ident = chosen.cluster))
    count.matrix <- as.data.frame(GetAssayData(tmp.s.obj.lv, slot = "data"))
    count.matrix <- count.matrix[gene_of_interest[[group]], ]
    countdf <- as.data.frame(rowMeans(count.matrix))
    colnames(countdf) <- c(chosen.cluster)
    countdf <- countdf %>% rownames_to_column("Gene")
    lv.df <- merge(lv.df, countdf, by.x = "Gene", by.y = "Gene")
  }
  lv.df <- lv.df %>% column_to_rownames("Gene")
  lv.df$total <- rowSums(lv.df)
  lv.df <- subset(lv.df, select = -c(total))
  
  pv.df <- data.frame(gene_of_interest[[group]])
  colnames(pv.df) <- c("Gene")
  for (chosen.cluster in cluster_of_interest){
    tmp.s.obj.pv <- subset(s.obj.PV, cells = WhichCells(object = s.obj.PV, ident = chosen.cluster))
    count.matrix <- as.data.frame(GetAssayData(tmp.s.obj.pv, slot = "data"))
    count.matrix <- count.matrix[gene_of_interest[[group]], ]
    countdf <- as.data.frame(rowMeans(count.matrix))
    colnames(countdf) <- c(chosen.cluster)
    countdf <- countdf %>% rownames_to_column("Gene")
    pv.df <- merge(pv.df, countdf, by.x = "Gene", by.y = "Gene")
  }
  pv.df <- pv.df %>% column_to_rownames("Gene")
  pv.df$total <- rowSums(pv.df)
  pv.df <- subset(pv.df, select = -c(total))
  
  colnames(lv.df) <- to_vec(for (item in colnames(lv.df)) sprintf("LV_%s", item))
  colnames(pv.df) <- to_vec(for (item in colnames(pv.df)) sprintf("PV_%s", item))
  
  lv.df <- lv.df %>% rownames_to_column("Gene")
  pv.df <- pv.df %>% rownames_to_column("Gene")
  
  merge.lv.pv.df <- merge(pv.df, lv.df, by.x = "Gene", by.y = "Gene") %>%
    column_to_rownames("Gene")
  
  merge.lv.pv.df$total <- rowSums(merge.lv.pv.df)
  merge.lv.pv.df <- subset(merge.lv.pv.df, merge.lv.pv.df$total != 0)
  merge.lv.pv.df <- subset(merge.lv.pv.df, select = -c(total))
  
  input.to.heatmap <- (merge.lv.pv.df - t(rowMeans(merge.lv.pv.df)))/(rowSds(as.matrix(merge.lv.pv.df)))
  tmp.heatmap.lv.pv <- input.to.heatmap %>% 
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-c(Gene), names_to = "Cluster", values_to = "Expression") %>%
    mutate(Cluster= fct_relevel(Cluster,colnames(input.to.heatmap))) %>%
    ggplot(aes(x=Cluster, y=Gene, fill=Expression)) + 
    geom_raster() + 
    scale_fill_gradient2(low="#1417f5", high="#f52314", guide="colorbar") +
    theme(text = element_text(size=25), axis.text.x = element_text(angle = 90))
  ggsave(plot = tmp.heatmap.lv.pv, filename = "HEATMAP_LV_PV.svg", path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)
  
  to.plot.lv.samples <- to_vec(for (item in colnames(input.to.heatmap)) if (grepl("LV", item)) item)
  to.plot.pv.samples <- to_vec(for (item in colnames(input.to.heatmap)) if (grepl("PV", item)) item)
  
  input.to.heatmap.lv.pv.only.lv <- input.to.heatmap[, to.plot.lv.samples]
  input.to.heatmap.lv.pv.only.pv <- input.to.heatmap[, to.plot.pv.samples]
  
  tmp.heatmap.lv <- input.to.heatmap.lv.pv.only.lv %>% 
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-c(Gene), names_to = "Cluster", values_to = "Expression") %>%
    mutate(Cluster= fct_relevel(Cluster,colnames(input.to.heatmap.lv.pv.only.lv))) %>%
    ggplot(aes(x=Cluster, y=Gene, fill=Expression)) + 
    geom_raster() + 
    scale_fill_gradient2(low="#1417f5", high="#f52314", guide="colorbar") +
    theme(text = element_text(size=25), axis.text.x = element_text(angle = 90))
  ggsave(plot = tmp.heatmap.lv, filename = "HEATMAP_LV.cross_condition.svg", path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)
  
  tmp.heatmap.pv <- input.to.heatmap.lv.pv.only.pv %>% 
    as.data.frame() %>%
    rownames_to_column("Gene") %>%
    pivot_longer(-c(Gene), names_to = "Cluster", values_to = "Expression") %>%
    mutate(Cluster= fct_relevel(Cluster,colnames(input.to.heatmap.lv.pv.only.pv))) %>%
    ggplot(aes(x=Cluster, y=Gene, fill=Expression)) + 
    geom_raster() + 
    scale_fill_gradient2(low="#1417f5", high="#f52314", guide="colorbar") +
    theme(text = element_text(size=25), axis.text.x = element_text(angle = 90))
  ggsave(plot = tmp.heatmap.pv, filename = "HEATMAP_PV.cross_condition.svg", path = path.to.save.output, device = "svg", width = 14, height = 10, dpi = 300)
  cat("\n \n")
}
