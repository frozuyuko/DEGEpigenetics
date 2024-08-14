library(maftools)
library(TCGAbiolinks)
library(dplyr)
tar_gz_file <- "TrainEpigenetics/TCGA-Assembler/gdc_download_20230529_140400.047300.tar.gz"
output_dir <- "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/ExtractedFiles"
untar(tar_gz_file, exdir = output_dir)
untar(tar_gz_file, exdir = "maf_file/TCGA-Assembler")
maf_mutect2 <- GDCquery_Maf("LUAD", pipelines = "mutect2") %>% subset(Sequencer == "Illumina HiSeq 2000") %>% read.maf
maf <- read.maf()