current_dir=$(pwd);

path_to_cellranger="/home/uk104163/cellranger/cellranger-7.2.0";
export PATH=${path_to_cellranger}:$PATH;

path_to_save_output=$HPCWORK;
mkdir -p ${path_to_save_output};
#path_to_ref_genome="/home/uk104163/references/refdata-gex-mm10-2020-A"
path_to_ref_genome="/home/uk104163/references/build-GRCh38/GRCh38";

path_to_fastq_files="/hpcwork/uk104163/raw_data/231205_Beckers_Berres_MedIII_fastq/compressed_tars/data/fastq/231205_A01742_0179_AHNWLNDRX3/mkfastq/outs/fastq_path/HNWLNDRX3";

for sample_id in Sample_3 Sample_4;do \
echo -e "Working on sample " $sample_id "\n";
cellranger multi --id=${sample_id} --csv=multi_config_${sample_id}.csv --localcores=12;done
