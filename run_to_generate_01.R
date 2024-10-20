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
path.to.storage <- "/media/hieunguyen/HD01/storage"
outdir <- "/home/hieunguyen/CRC1382/outdir"
path.to.rmd <- file.path(path.to.project.src, "01_preliminary_analysis.Rmd")

for (batch.id in c("batch1", "batch2")){
  save.html.name <- str_replace(basename(path.to.rmd), ".Rmd", sprintf(".%s.html", batch.id))
  path.to.save.html <- file.path(outdir, PROJECT, batch.id, "html_outputs")
  dir.create(path.to.save.html, showWarnings = FALSE, recursive = TRUE)
  
  rmarkdown::render(input = path.to.rmd, 
                    params = list(
                      batch.id = batch.id
                    ),
                    output_file = save.html.name, 
                    output_dir = path.to.save.html)
}

