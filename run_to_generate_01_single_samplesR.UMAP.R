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
PROJECT <- "ABeckers_data_single_samples"
path.to.storage <- "/media/hieunguyen/HD01/storage"
outdir <- "/home/hieunguyen/CRC1382/outdir"
path.to.rmd <- file.path(path.to.project.src, "01_preliminary_analysis.single_sample_UMAP.Rmd")

for (output.version in c("v0.1", "default")){
  for (sample.id in c("LV20", "PV20", "LV21", "PV21", "LV73", "PV73")){
    save.html.name <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", sample.id))
    path.to.save.html <- file.path(outdir, PROJECT, output.version, sample.id, "html_outputs")
    dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
    
    if (file.exists(file.path(path.to.save.html, save.html.name)) == FALSE){
      print(sprintf("Sample ID %s", sample.id))
      print(sprintf("output version: %s", output.version))
      
      rmarkdown::render(input = path.to.rmd, 
                        params = list(
                          sample.id = sample.id,
                          output.version = output.version
                        ),
                        output_file = save.html.name, 
                        output_dir = path.to.save.html)
    }
  }
}

