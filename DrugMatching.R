drug <- read.csv("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/DrugLIstLUAD.csv", header = TRUE)
View(drug)
#Gene Data
View(early.alive)
drug.match <- intersect(colnames(early.alive), drug$bcr_patient_barcode)
colnames(early.alive)
View(early.alive)
drug.match <- as.data.frame(drug.match)
View(drug.match)
colnames(drug.match)[colnames(drug.match) == "drug.match"] <- "bcr_patient_barcode"
View(drug.match)

barcoding <- intersect(drug$bcr_patient_barcode, drug.match$bcr_patient_barcode)

subset_data <- drug[drug$bcr_patient_barcode %in% barcoding, ]
subset_data

drugcombined <- merge(drug.match, subset_data[, c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name")], 
                      by = "bcr_patient_barcode", all.x = TRUE)
View(drugcombined)

drugcombined$pharmaceutical_therapy_drug_name
drugname <- unique(drugcombined$pharmaceutical_therapy_drug_name)
drugname <- as.data.frame(drugname)
View(drugname)

#miRNA
drugmiRNA <- read.csv("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/top11updownregulated.csv")
drugmiRNA <- drugmiRNA[drugmiRNA$status == "cancer", ]
View(drugmiRNA)
drugmiRNA.match <- intersect(drugmiRNA$Barcode, drug$bcr_patient_barcode)
View(drugmiRNA.match)
drugmiRNA.match <- as.data.frame(drugmiRNA.match)
colnames(drugmiRNA.match)[colnames(drugmiRNA.match) == "drugmiRNA.match"] <- "bcr_patient_barcode"
drugmiRNA.match
miRNAbarcoding <- intersect(drug$bcr_patient_barcode, drugmiRNA.match$bcr_patient_barcode)
miRNAbarcoding
datasubset <- drug[drug$bcr_patient_barcode %in% miRNAbarcoding, ]
combineddrug <- merge(drugmiRNA.match, datasubset[, c("bcr_patient_barcode", "pharmaceutical_therapy_drug_name")], 
                      by = "bcr_patient_barcode", all.x = TRUE)
View(combineddrug)
combineddrug$pharmaceutical_therapy_drug_name
namedrug <- unique(combineddrug$pharmaceutical_therapy_drug_name)
namedrug <- as.data.frame(namedrug)
View(namedrug)
View(drugname)
colnames(namedrug)[colnames(namedrug) == "namedrug"] <- "drugname"
library(dplyr)
uniquedrug <- anti_join(namedrug, drugname, by = "drugname")
print(uniquedrug)
print(unique_df1)
drugmiRNAgene <- append(drugname$drugname, namedrug$drugname)
drugmiRNAgene <- unique(drugmiRNAgene)
drugmiRNAgene <- as.data.frame(drugmiRNAgene)
View(drugmiRNAgene)
drugmiRNAgene

#Intersect barcode
druganalysis <- unique(drugcombined$bcr_patient_barcode, combinedrug)
combinedrug

count <- sum(ifelse(grepl("Paclitaxel|Taxol", drugcombined$pharmaceutical_therapy_drug_name, 1,0)))