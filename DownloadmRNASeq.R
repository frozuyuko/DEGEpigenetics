setwd("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler")
getwd()
RNAExp <- DownloadRNASeqData(cancerType = "LUAD", assayPlatform =
                                   "gene.normalized_RNAseq",
                                 saveFolderName = "DownloadRNASeqData")
list_RNAExp <- ProcessRNASeqData(inputFilePath =
                                       RNAExp[1],outputFileName = "ProcessRNASeqData", 
                                     outputFileFolder =
                                       "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/",
                                 dataType	=	"geneExp",	verType	=	"RNASeqV2")
listRNAExp <- ProcessRNASeqData(inputFilePath =
                                   RNAExp[1],outputFileName = "ProcessRNASeqDataTest", 
                                 outputFileFolder =
                                   "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/",
                                 dataType	=	"geneExp")

#loading gene data
det <- load("ProcessRNASeqData.rda")
View(head(det))
View(head(Data))
Des[,1] <- apply(Des[,1:2], 1, function(x) paste(x, collapse="|"))
View(head(Des))
colnames(Des)[1] <-"gene_id"
gene_id <- Des
row.names(Data) <- Des
View(head(Data))
temp <- Data
View(temp)
#trash
RNAexpress <- read.csv(file = "RNAseq.csv")
View(RNAexpress)
GeneName <- read.table(file = "GeneSymbol.csv")
colnames(GeneName)[1] <- "gene_name"
View(GeneName)
rownames(GeneName) <- RNAexpress
RNAexpress 
patient.type <- sapply(colnames(RNAexpress), function(s) unlist(strsplit(s, "-"))[4])
tumor.primary <- RNAexpress[, grep("^01", patient.type)]
normale <- RNAexpress[, c(grep("^10", patient.type), grep("^11", patient.type), grep("^12", patient.type))]
tum.luad <- substr(colnames(tumor.primary), 1, 12) 
View(tum.luad)
norm.luad <- substr(colnames(normale), 1, 12)
colnames(tumor.primary) <- tum.luad
View(head(tumor.primary))

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
BiocManager::install("TCGAassembler")
if (!require("BiocManager", quietly = TRUE))
  install.packages("TCGAassembler")

query <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Transcriptome Profiling", 
                  data.type = "Gene Expression Quantification",
                  legacy = TRUE)
GDCdownload(query)
