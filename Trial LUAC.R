# to browse: https://portal.gdc.cancer.gov/ 
# about tcga: https://github.com/compgenome365/TCGA-Assembler-2
# to understand assayplatform: https://github.com/compgenome365/TCGA-Assembler-2/blob/master/TCGA-Assembler/Module_A.R
setwd("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler")
getwd()
install.packages("httr")
install.packages("HGNChelper")
install.packages("RCurl")
install.packages("rjson")
install.packages("stringr")
library(httr)
library(HGNChelper)
library(RCurl)
library(rjson)
library(stringr)
source("./Module_A.R")
source("./Module_B.R")
miRNAExp <- DownloadmiRNASeqData(cancerType = "LUSC", assayPlatform =
                                   "mir_HiSeq.hg19.mirbase20",
                                 saveFolderName = "DownloadmiRNASeqData")
miRNAIsoformExpLung <- DownloadmiRNASeqData(cancerType = "LUSC", assayPlatform =
                                       "mirIsoform_HiSeq.hg19.mirbase20",
                                       saveFolderName = "DownloadmiRNAisoformSeqLungData")
list_miRNAExp <- ProcessmiRNASeqData(inputFilePath =
                                       miRNAExp[1],outputFileName = "ProcessmiRNASeqData", 
                                     outputFileFolder =
                                       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/")
LUAD_clinical <- DownloadBiospecimenClinicalData(cancerType = "LUAD",
                                                 saveFolderName = 
                                                   "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/DownloadBiospecimenClinicalData")

dat <- load("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/ProcessmiRNASeqData__RPM.rda")
View(head(dat))
View(head(Data))
row.names(Data) <- Des
View(head(Des))
View(head(Data))
temp <- Data
View(temp)
temp.type <- sapply(colnames(temp), function(s) unlist(strsplit(s, "-"))[4])
primary.tumor <- temp[, grep("^01", temp.type)]
normal <- temp[, c(grep("^10", temp.type), grep("^11", temp.type), grep("^12", temp.type))]
tum.bcr <- substr(colnames(primary.tumor), 1, 12) 
norm.bcr <- substr(colnames(normal), 1, 12)
colnames(primary.tumor) <- tum.bcr
View(head(primary.tumor))
colnames(normal) <- norm.bcr
View(head(normal))
install.packages("dplyr")
install.packages("tidyverse")
library(dplyr)
library(tidyverse)
library("readxl")
clinical <- read_excel("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/DownloadBiospecimenClinicalData/ClinicalPatientLUAD.xlsx")
clinical <- read.csv("ClinicalPatientLUAD.csv")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
if (!require("edgeR", quietly = TRUE))
  install.packages("edgeR")

BiocManager::install("HTSFilter")
install.packages("limma")
install.packages("edgeR")

install.packages("HTSFilter")
library(limma)
library(edgeR)
library(HTSFilter)
rawCountTable = load("blackS1_data.rda") 

RNAExp <- DownloadRNASeqData(cancerType = "LUAD", assayPlatform =
                                   "gene.normalized_RNAseq",saveFolderName = "DownloadRNASeqData")
list_RNAExp <- ProcessRNASeqData(inputFilePath =
                                       RNAExp[1],outputFileName = "ProcessRNASeqData", 
                                     outputFileFolder =
                                       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/")

dat <- load("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/ProcessRNASeqData__RPM.rda")
View(head(Data))
row.names(Data) <- Des
View(head(Data))
temp <- Data
DNAExp <- DownloadMethylationData(cancerType = "LUAD", assayPlatform =
                               "methylation_450",saveFolderName = "DownloadMethylationData")

#Early cancer data
stage_1.barcode <- clinical[(clinical$ajcc_pathologic_tumor_stage == "Stage I" |
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA" | clinical$ajcc_pathologic_tumor_stage == "Stage IA1"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA2"| clinical$ajcc_pathologic_tumor_stage == "Stage IB"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IB1"| clinical$ajcc_pathologic_tumor_stage == "Stage IB2"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IC"| clinical$ajcc_pathologic_tumor_stage == "Ia" |
                               clinical$ajcc_pathologic_tumor_stage == "Ib"| clinical$ajcc_pathologic_tumor_stage == "I"), ]
stage_2.barcode <- clinical[(clinical$ajcc_pathologic_tumor_stage == "Stage II"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IIA" | clinical$ajcc_pathologic_tumor_stage == "Stage IIA1"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IIA2" | clinical$ajcc_pathologic_tumor_stage == "Stage IIC"|
                               clinical$ajcc_pathologic_tumor_stage == "IIa" | clinical$ajcc_pathologic_tumor_stage == "IIb"|
                               clinical$ajcc_pathologic_tumor_stage == "II"), ]
stage_12.barcode <- clinical[(clinical$ajcc_pathologic_tumor_stage == "Stage I" |
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA" | clinical$ajcc_pathologic_tumor_stage == "Stage IA1"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA2"| clinical$ajcc_pathologic_tumor_stage == "Stage IB"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IB1"| clinical$ajcc_pathologic_tumor_stage == "Stage IB2"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IC"| clinical$ajcc_pathologic_tumor_stage == "Ia" |
                               clinical$ajcc_pathologic_tumor_stage == "Ib"| clinical$ajcc_pathologic_tumor_stage == "I" clinical$ajcc_pathologic_tumor_stage == "Stage I" |
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA" | clinical$ajcc_pathologic_tumor_stage == "Stage IA1"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IA2"| clinical$ajcc_pathologic_tumor_stage == "Stage IB"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IB1"| clinical$ajcc_pathologic_tumor_stage == "Stage IB2"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IC"| clinical$ajcc_pathologic_tumor_stage == "Ia" |
                               clinical$ajcc_pathologic_tumor_stage == "Ib"| clinical$ajcc_pathologic_tumor_stage == "I"| 
                               clinical$ajcc_pathologic_tumor_stage == "Stage II"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IIA" | clinical$ajcc_pathologic_tumor_stage == "Stage IIA1"|
                               clinical$ajcc_pathologic_tumor_stage == "Stage IIA2" | clinical$ajcc_pathologic_tumor_stage == "Stage IIC"|
                               clinical$ajcc_pathologic_tumor_stage == "IIa" | clinical$ajcc_pathologic_tumor_stage == "IIb"|
                               clinical$ajcc_pathologic_tumor_stage == "II"), ]
#Alive data
alive.barcode <- clinical[(clinical$vital_status == "Alive"), ]
alive.digit <- alive.barcode$bcr_patient_barcode

stage_1.barcode <- stage_1.barcode$bcr_patient_barcode
stage_2.barcode <- stage_2.barcode$bcr_patient_barcode

#normal data
normal_alive <- intersect(norm.bcr, alive.digit)
normal_alive
dat_normal <- normal [, normal_alive]
View(dat_normal)
write.csv(dat_normal, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/normalalive.csv")
save(dat_normal, file = 
       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/normalalive.rda")

#cancer data
cancer_alive <- intersect(tum.bcr, alive.digit)
dat_cancer <- primary.tumor [, cancer_alive]
write.csv(dat_cancer, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/canceralive.csv")
save(dat_cancer, file = 
       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/canceralive.rda")

#specific stage
stage_1.alive <- intersect(cancer_alive, stage_1.barcode)
stage_2.alive <- intersect(cancer_alive, stage_2.barcode)
stage_1.match <- primary.tumor [, stage_1.alive]
stage_2.match <- primary.tumor [, stage_2.alive]

clinical_early.tumor <- intersect()

early.alive <- cbind(stage_1.match, stage_2.match)
early.match <- intersect(colnames(early.alive), clinical$bcr_patient_barcode)
early.match
earlytumor.match <- intersect(tum.bcr, early.match)
earlytumor.match
earlytumor.match_norm <- intersect(norm.bcr, colnames(early.match))
#upregulate and downregulate
rawCount <- cbind(primary.tumor[, colnames(early.alive)],
                  normal[, colnames(dat_normal)])
rawCount <- rawCount
rawCountresult <- unique(rawCount, by = colnames(rawCount))
colnames(rawCountresult) <- make.unique(colnames(rawCountresult), sep = "_")
rawCountresult <- t(rawCountresult)
View(rawCountresult)
earlycancer_alive <- t(rawCountresult)
View(earlycancer_alive)
write.csv(earlycancer_alive, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/miRNALUAD.csv")

early.match <- as.data.frame(early.match)
View(early.match)
early.match$tumor <- "cancer"
normal_alive <- as.data.frame(normal_alive)
normal_alive$tumor <- "normal"
View(normal_alive)
early.match
normal_alive
colnames(early.match)[1] <- "patient"
sampleInfo <- rbind(early.match, normal_alive)
View(sampleInfo)
colnames(sampleInfo)[2] <- "status"

# create DGEList object
dgeFull <- DGEList(rawCountresult, group=sampleInfo$tumor)
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
dgeFull <- calcNormFactors(dgeFull, method="TMM")
dgeFull$samples
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

dgeFull
View(dgeFull)
dgeTest <- exactTest(dgeFull)
dgeTest <- t(dgeTest)
dgeTest

dgeframe <- as.data.frame(dgeTest)

install.packages("HTSFilter")
library(HTSFilter)

par(mar=c(5, 4, 4, 2) + 0.1)
#adjust plot margins
par(mar = c(1, 1, 1, 1))

#create scatterplot
plot(1:30)
filtData <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtData)
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table), sort.by = "logFC", p.value = 0.01)

head(resFilt)

sigDownReg <- resFilt$table

sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[sigDownReg$logFC<(-1.5),]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
sigDownReg <- sigDownReg[1:10, ]

sigUpReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigUpReg <- sigUpReg[sigUpReg$logFC>1.5,]
sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing = TRUE),]
sigUpReg <- sigUpReg[1:10, ]

sigDEG <- rbind(sigDownReg, sigUpReg)

#Problem hereee::: 
alive_tumor <- primary.tumor[, early.match]
alive_tumor <- t(alive_tumor)
down.express <- alive_tumor[rownames(sigDownReg),]
up.express <- alive_tumor[rownames(sigUpReg),]
ab.express <- alive_tumor[rownames(sigDEG),]

down.express <- down.express[apply(down.express == 0, 1, sum) <= 150,]
up.express <- up.express[apply(up.express == 0, 1, sum) <= 150,]
ab.express <- ab.express[apply(ab.express == 0, 1, sum) <= 150,]

cur_res <- decideTestsDGE(dgeTest, adjust.method = "BH", p.value = 0.05)
cur_res

sel_deg <- which(cur_res[ ,1] != 0)
View(cur_res)
length(sel_deg)

DEG_pairwiseExact <- rownames(cur_res)[sel_deg]
head(DEG_pairwiseExact)
DEG_pairwiseExact


write.csv(stage_1.match, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/cancerstage1alive.csv")
save(stage_1.match, file = 
       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/cancerstage1alive.rda")
write.csv(stage_2.match, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/cancerstage2alive.csv")
save(stage_2.match, file = 
       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/cancerstage2alive.rda")


cancer_combined <- primary.tumor [, stage_1.alive] + primary.tumor [, stage_2.alive]
cancer_combined <- bind_rows(stage_1.match, stage_2.match)

#packages
library(TCGA2STAT)
library(dplyr)
stage1_filtered <- filter(stage_1.match, genes %in% intersection_data)
stage2_filtered <- filter(stage_2.match, genes %in% intersection_data)
cancer_combined <- bind_rows(stage_1.match, stage_2.match)

write.csv(dat_cancer, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/canceralive.csv")
save(dat_cancer, file = 
       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/canceralive.rda")

#convert cancer data


cancer.match <- intersect(stage_1.match, stage_2.match)
cancernormal.match <- intersect(cancer.match, normal_alive)







# Ethnic:
indian.barcode <- clinical[(clinical$race == "AMERICAN INDIAN OR ALASKA NATIVE"), ]
asian.barcode <- clinical[(clinical$race == "ASIAN"), ]
black.barcode <- clinical[(clinical$race == "BLACK OR AFRICAN AMERICAN"), ]
white.barcode <- clinical[(clinical$race == "WHITE"), ]
white_race.barcode <- clinical[(clinical$race == "WHITE" |
                                  clinical$race == "white" | clinical$race == "White"), ]
white.digit <- white.barcode$bcr_patient_barcode
white.match_tum<- intersect(tum.bcr, white.digit)
white.match_nor<- intersect(norm.bcr, white.digit)


all.barcode <- clinical[(clinical$race == "WHITE"),
                        (clinical$race == "AMERICAN INDIAN OR ALASKA NATIVE"),
                        (clinical$race == "ASIAN"),]

stage_1.match <- intersect(tum.bcr, stage_1.barcode)
stage_2.match <- intersect(tum.bcr, stage_2.barcode)

stage_1.alive <- intersect(stage_1.match, alive.barcode)
stage_2.alive <- intersect(tum.bcr, stage_2.barcode)

stage1.dat <- primary.tumor[, stage_1.match]
stage2.dat <- primary.tumor[, stage_2.match]


white.dat_tumor <- primary.tumor [, white.match_tum]
white.dat_tumor_stage1_barcode <- intersect (stage_1.match, white.match_tum)
white.dat_tumor_stage1_dataset <- primary.tumor [, white.dat_tumor_stage1_barcode]

write.csv(whiteTumor.dat_tumor, file = "D:\\Training Cancer Epigenetic\\TCGAAssembler2\\whiteTumordata.csv")
#e.g. save whiteTumor.dat_tumor in .rda
save(whiteTumor.dat_tumor, file = "D:\\Training Cancer Epigenetic\\TCGAAssembler2\\whiteTumordata.rda")

white_race.match_norm <- intersect(norm.bcr, white.barcode)






---------

#Support Vector Machines
install.packages('caTools')
library(caTools)
install.packages('e1071')
library(e1071)
clinicaldataset = clinical[3:5]
dataset = dataset[3:5]
clinicaldata = clinical[]
install.packages('Chemmine')
library()