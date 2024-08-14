det <- load("ProcessRNASeqDataTest.rda")
View(head(det))
View(head(Data))
View(head(Des))
# Replace "filename.csv" with the actual path and name of your CSV file
GeneList <- read.csv("EntrezIDRNASeq.csv",header = TRUE)
GeneList <- GeneList[1]
View(GeneList)
Data <- t(Data)
View(Data)
colnames(GeneList) <- Data
colnames(Data) <- GeneList$GeneSymbol
View(GeneList)
nrow(Data)
nrow(GeneList)
length(Data)
as.vector(GeneList)
row.names(GeneList) <- Data
GeneList$GeneSymbol
rownames(GeneList$GeneSymbol) <- Data
rownames(GeneList) <- 
rownames(Data) <- GeneList[,1]
View(Data)
row.names(Data) <- GeneList
Des <- colnames(GeneList)[2]
View(Des)
Gene <- t(Data)
View(Gene)
Gene <- t(Gene)
patient.type <- sapply(colnames(Gene), function(s) unlist(strsplit(s, "-"))[4])
tumor.primary <- Gene[, grep("^01", patient.type)]
normal <- Gene[, c(grep("^10", patient.type), grep("^11", patient.type), grep("^12", patient.type))]
tum.luad <- substr(colnames(tumor.primary), 1, 12) 
norm.luad <- substr(colnames(normal), 1, 12)
colnames(tumor.primary) <- tum.luad
View(head(tumor.primary))
colnames(normal) <- norm.luad
View(head(normal))

#normal data
alivenormal <- intersect(norm.luad, alive.digit)
normaldat <- normal [, alivenormal]
View(normaldat)

#cancer data
cancer_alive <- intersect(tum.luad, alive.digit)
tum.luad
dat_cancer <- tumor.primary[, cancer_alive]

#specific stage
stage_1.alive <- intersect(cancer_alive, stage_1.barcode)
stage_2.alive <- intersect(cancer_alive, stage_2.barcode)
stage_1.match <- tumor.primary [, stage_1.alive]
stage_2.match <- tumor.primary [, stage_2.alive]

early.alive <- cbind(stage_1.match, stage_2.match)
earlymatch <- intersect(colnames(early.alive), clinical$bcr_patient_barcode)
earlymatch
earlytumor.match <- intersect(tum.luad, earlymatch)
earlytumor.match_norm <- intersect(norm.luad, colnames(normaldat))
View(earlytumor.match_norm)
#upregulate and downregulate
rawCount <- cbind(tumor.primary[, colnames(early.alive)],
                  normal[, colnames(normaldat)])
View(rawCount)
rawCountresult <- unique(rawCount, by = colnames(rawCount))
colnames(rawCountresult) <- make.unique(colnames(rawCountresult), sep = "_")
View(rawCountresult)
rawCountresult <- t(rawCountresult)
earlymatch <- as.data.frame(earlymatch)
earlymatch$Condition <- "cancer"
alivenormal <- as.data.frame(alivenormal)
alivenormal$Condition <- "normal"
colnames(alivenormal)[1] <- "patient"
colnames(earlymatch)[1] <- "patient"
View(alivenormal)
sampleinformation <- rbind(earlymatch, alivenormal)

View(sampleinformation)
write.csv(sampleInfo, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/sampleInfo.csv")
# create DGEList object
dgeFull <- DGEList(rawCountresult, group=sampleInfo$Condition)
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
                   group=dgeFull$samples$group)
dgeFull <- calcNormFactors(dgeFull, method="RLE")
dgeFull$samples
dgeFull <- estimateCommonDisp(dgeFull)
dgeFull <- estimateTagwiseDisp(dgeFull)

dgeTest <- exactTest(dgeFull)
library(HTSFilter)
par(mar=c(5, 4, 4, 2) + 0.1)
#adjust plot margins
par(mar = c(1, 1, 1, 1))

#create scatterplot
plot(1:30)
filtRLE <- HTSFilter(dgeFull)$filteredData
dgeTestFilt <- exactTest(filtRLE)
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTestFilt$table), sort.by = "logFC", p.value = 0.01)
head(resFilt)

sigDownReg <- resFilt$table

sigDownReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigDownReg <- sigDownReg[sigDownReg$logFC<(-2.5),]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
sigDownReg <- sigDownReg[1:10, ]
sigDownReg

sigUpReg <- resFilt$table[resFilt$table$FDR<0.01,]
sigUpReg <- sigUpReg[sigUpReg$logFC>1.5,]
sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing = TRUE),]
sigUpReg <- sigUpReg[1:10, ]
sigUpReg

sigUpDown <- rbind(sigDownReg, sigUpReg)
sigUpDown

#Problem hereee::: 
alive_tumor <- tumor.primary[, earlytumor.match]
down.express <- alive_tumor[rownames(sigDownReg),]
up.express <- alive_tumor[rownames(sigUpReg),]
updown.express <- alive_tumor[rownames(sigUpDown),]

down.express <- down.express[apply(down.express == 0, 1, sum) <= 82,]
View(down.express)
up.express <- up.express[apply(up.express == 0, 1, sum) <= 82,]
View(up.express)
updown.express <- updown.express[apply(updown.express == 0, 1, sum) <= 82,]
View(updown.express)
updown.express
up.express <- as.data.frame(t(up.express))
View(up.express)
up.express$status_up <- "upregulated"
down.express <- as.data.frame(t(down.express))
View(down.express)
down.express$status_down <- "downregulated"
up.express$status <- "cancer"
down.express$status <- "cancer"
updown.express <- as.data.frame(t(updown.express))
updown.express$status <- "cancer"
View(updown.express)

expression <- bind_cols(up.express, down.express)
expression$status
View(expression)
colnames(expression)[11] <- "status_up"
colnames(expression)[21] <- "status_down"
View(expression)


write.csv(combinedupdown, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/combinedupdown.csv")

#extract normal data
# Identify the column names
normal.data <- t(dat_normal)
View(normal.data)
normalalive <- colnames(normal.data)
normalalive
updowncombined <- colnames(expression)
updowncombined
updowncomb <- colnames(updown.express)
upcombined <- colnames(up.express)
downcombined <- colnames(down.express)

# Extract the columns
normalupdown <- intersect(normalalive, updowncombined)
View(normalupdown)
normalup <- intersect(normalalive, upcombined)
normaldown <- intersect(normalalive, downcombined)
updownnormal <- intersect(normalalive, updowncomb)
View(updownnormal)

#recall the data
normal_intersect <- normal.data[, normalupdown] 
normal_intersect <- as.data.frame(normal_intersect)
normal_upintersect <- as.data.frame(normal.data[, normalup])
normal_downintersect <- as.data.frame(normal.data[, normaldown])
normal_updownintersect <- as.data.frame(normal.data[, updownnormal])
View(normal_updownintersect)
position <- 15
# Add an empty column at the specified position
normal_intersect <- cbind(normal_intersect[, 1:(position-1)], "", 
                          normal_intersect[, position:ncol(normal_intersect)])
normal_intersect$status_up <- "normal"
normal_updownintersect$status <- "normal"
normal_upintersect$status <- "normal"
normal_downintersect$status <- "normal"
View(normal_updownintersect)

# Rename the new column if needed
colnames(normal_intersect)[position] <- "status_up"

# View the updated data frame
print(normal_upintersect)

normal_intersect[, 11] <- "normal"


View(normal_intersect)


#combined upregulate and normal
View(combinedupdown)
updownregulate_normal <- rbind(normal_intersect, expression)
upregulate_normal <- rbind(normal_upintersect, up.express)
View(upregulate_normal)
downregulate_normal <- rbind(normal_downintersect, down.express)
View(downregulate_normal)
normal_updownregulate <- rbind(normal_updownintersect, updown.express)
View(normal_updownregulate)

#save as csv
write.csv(normal_updownregulate, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/top14updowngenes.csv")
write.csv(upregulate_normal, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/top9upgenes.csv")
write.csv(downregulate_normal, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/top5downgenes.csv")
