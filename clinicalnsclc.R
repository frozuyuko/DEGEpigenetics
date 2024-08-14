#aggregate clinical data 
clinicalluad <- read.csv("ClinicalPatientLUAD.csv")
clinicallusc <- read.csv("LUSCclinical.csv", sep=";")
clinicalluadtcga <- clinicalluad[c("bcr_patient_barcode","ajcc_pathologic_tumor_stage", "vital_status")]
colnames(clinicalluadtcga)[1] = "barcode"
colnames(clinicalluadtcga)[2] = "pathologic_stage"
clinicallusctcga <- clinicallusc[c("submitter_id","ajcc_pathologic_stage", "vital_status")]
colnames(clinicallusctcga)[1] = "barcode"
colnames(clinicallusctcga)[2] = "pathologic_stage"
clinicalnsclc <- rbind(clinicalluadtcga,clinicallusctcga)
View(clinicalnsclc)

#Obtain early stages (stage I and II)
# Early cancer data
stage_12.barcode <- clinicalnsclc[(clinicalnsclc$pathologic_stage == "Stage I" |
                                     clinicalnsclc$pathologic_stage == "Stage IA" | clinicalnsclc$pathologic_stage == "Stage IA1"|
                                     clinicalnsclc$pathologic_stage == "Stage IA2"| clinicalnsclc$pathologic_stage == "Stage IB"|
                                     clinicalnsclc$pathologic_stage == "Stage IB1"| clinicalnsclc$pathologic_stage == "Stage IB2"|
                                     clinicalnsclc$pathologic_stage == "Stage IC"| clinicalnsclc$pathologic_stage == "Ia" |
                                     clinicalnsclc$pathologic_stage == "Ib"| clinicalnsclc$pathologic_stage == "I" | clinicalnsclc$pathologic_stage == "Stage I" |
                                     clinicalnsclc$pathologic_stage == "Stage IA" | clinicalnsclc$pathologic_stage == "Stage IA1"|
                                     clinicalnsclc$pathologic_stage == "Stage IA2"| clinicalnsclc$pathologic_stage == "Stage IB"|
                                     clinicalnsclc$pathologic_stage == "Stage IB1"| clinicalnsclc$pathologic_stage == "Stage IB2"|
                                     clinicalnsclc$pathologic_stage == "Stage IC"| clinicalnsclc$pathologic_stage == "Ia" |
                                     clinicalnsclc$pathologic_stage == "Ib"| clinicalnsclc$pathologic_stage == "I"| 
                                     clinicalnsclc$pathologic_stage == "Stage II"|
                                     clinicalnsclc$pathologic_stage == "Stage IIA" | clinicalnsclc$pathologic_stage == "Stage IIA1"|
                                     clinicalnsclc$pathologic_stage == "Stage IIA2" | clinicalnsclc$pathologic_stage == "Stage IIC"|
                                     clinicalnsclc$pathologic_stage == "IIa" | clinicalnsclc$pathologic_stage == "IIb"|
                                     clinicalnsclc$pathologic_stage == "II"), ]
# Alive data
alive.barcode <- clinicalnsclc[(clinicalnsclc$vital_status == "Alive"), ]
alive.digit <- alive.barcode$barcode
stage_12.barcode <- stage_12.barcode$barcode

#miRNA LUAD 

#miRNA LUSC 


#mRNA LUAD 
geneluad <- read.table("ProcessRNASeqDataTest.txt", header=TRUE)
geneluad$gene_id <- paste(geneluad$`GeneSymbol`, geneluad$`EntrezID`, sep = "|")
geneluad <- geneluad[, !(names(geneluad) %in% c("GeneSymbol", "EntrezID"))]
geneluad <- t(geneluad)
colnames(geneluad) <- geneluad[577,]
geneluad <- geneluad[!(rownames(geneluad) == "gene_id"), ]
#mRNA LUSC
genelusc <- read.table("LUSC__gene.normalized_RNAseq__tissueTypeAll__20230509133651.txt", header=TRUE)
rownames(genelusc) <- genelusc[, 1]
genelusc <- genelusc[, -1]
genelusc <- t(genelusc)
#NSCLC mRNA
write.csv(geneluad, "geneluad.csv", row.names = TRUE)
write.csv(genelusc, "genelusc.csv", row.names = TRUE)
genensclc <- read.csv("genensclc.csv")
genensclc <- genensclc[,-1]
rownames(genensclc) <- genensclc[,1]
genensclc <- genensclc[,-1]

#miRNA LUAD 
mirnaluad <- read.table("ProcessmiRNASeqDAta__RPM.txt", sep="\t", header = TRUE)
rownames(mirnaluad) <- mirnaluad[,1]
mirnaluad <- mirnaluad[,-1]
#miRNA LUSC 
mirnalusc <- read.table()






# above is trash

#miRNA NSCLC
mirnansclc <- read.table('NSCLC__mir_HiSeq.hg19.mirbase20__tissueTypeAll.txt', header = TRUE)
mirnansclc <- mirnansclc[-1,]
rownames(mirnansclc) <- mirnansclc[,1]
mirnansclc <- mirnansclc[,-1]
#Repair mirna data 
mirnansclc.type <- sapply(colnames(mirnansclc), function(s) unlist(strsplit(s, "\\."))[4])
#mRNA NSCLC
genensclc <- read.table('NSCLC__gene.normalized_RNAseq__tissueTypeAll.txt', header = TRUE)
rownames(genensclc) <- genensclc[,1]
genensclc <- genensclc[,-1]
#Repair mRNA data 
genensclc.type <- sapply(colnames(genensclc), function(s) unlist(strsplit(s, "\\."))[4])
#match samples with mirna 
#overallnsclc.type <- intersect(colnames(genensclc),colnames(mirnansclc))
nsclcgene.primarytumor <- genensclc[, grep("^01", genensclc.type)]
nsclcgene.normal <- genensclc[, c(grep("^10", genensclc.type), grep("^11", genensclc.type), grep("^12", genensclc.type))]
genensclc.tumor <- substr(colnames(nsclcgene.primarytumor), 1, 12)
genensclc.normal <- substr(colnames(nsclcgene.normal), 1, 12)
genensclc.tumor <- gsub("\\.", "-", genensclc.tumor)
genensclc.normal <- gsub("\\.", "-", genensclc.normal)
#Drop the 5 if it appears 
#genensclc.tumor <- gsub("\\.\\d+$", "", genensclc.tumor)
#genensclc.normal <- gsub("\\.\\d+$", "", genensclc.normal) 
colnames(nsclcgene.primarytumor) <- genensclc.tumor
colnames(nsclcgene.normal) <- genensclc.normal
#do the same for gene
nsclcmirna.primarytumor <- mirnansclc[, grep("^01", mirnansclc.type)]
nsclcmirna.normal <- mirnansclc[, c(grep("^10", mirnansclc.type), grep("^11", mirnansclc.type), grep("^12", mirnansclc.type))]
mirnansclc.tumor <- substr(colnames(nsclcmirna.primarytumor), 1, 12)
mirnansclc.normal <- substr(colnames(nsclcmirna.normal), 1, 12)
mirnansclc.tumor <- gsub("\\.", "-", mirnansclc.tumor)
mirnansclc.normal <- gsub("\\.", "-", mirnansclc.normal)
#Drop the 5 if it appears 
#mirnansclc.tumor <- gsub("\\.\\d+$", "", mirnansclc.tumor)
#mirnansclc.normal <- gsub("\\.\\d+$", "", mirnansclc.normal) 
colnames(nsclcmirna.primarytumor) <- mirnansclc.tumor
colnames(nsclcmirna.normal) <- mirnansclc.normal
#matching barcode 
overalltumor <- intersect(colnames(nsclcgene.primarytumor), colnames(nsclcmirna.primarytumor))
nsclcgene.primarytumor <- nsclcgene.primarytumor[, overalltumor]
nsclcmirna.primarytumor <- nsclcmirna.primarytumor[, overalltumor]
overallnormal <- intersect(colnames(nsclcgene.normal), colnames(nsclcmirna.normal))
nsclcgene.normal <- nsclcgene.normal[, overallnormal]
nsclcmirna.normal <- nsclcmirna.normal[, overallnormal]

#normal alive 
overallnormalalive <- intersect(genensclc.normal, alive.digit)

#tumor gene alive
overalltumoralive <- intersect(genensclc.tumor, alive.digit)

#only early stage patients
overalltumorearlyalive <- intersect(overalltumoralive,stage_12.barcode)
overalltumorearlyalive <- intersect(overalltumorearlyalive, alive.digit)

#data preparation
#Early cancer
#tumorgenenames <- setNames(rep(TRUE, length(colnames(nsclcgene.primarytumor))), colnames(nsclcgene.primarytumor))
#genetumormatches <- overalltumorearlyalive %in% names(tumorgenenames)
#matching_genetumor <- overalltumorearlyalive[genetumormatches]
earlytumorgene.match <- nsclcgene.primarytumor[, matching_genetumor]
#earlytumorgene.match <- nsclcgene.primarytumor[, overalltumorearlyalive]
tumormirnanames <- setNames(rep(TRUE, length(colnames(nsclcmirna.primarytumor))), colnames(nsclcmirna.primarytumor))
mirnatumormatches <- overalltumorearlyalive %in% names(tumormirnanames)
matching_mirnatumor <- overalltumorearlyalive[mirnatumormatches]
earlytumormirna.match <- nsclcmirna.primarytumor[, matching_mirnatumor]
#Normal 
normalgenenames <- setNames(rep(TRUE, length(colnames(nsclcgene.normal))), colnames(nsclcgene.normal))
genenormalmatches <- overallnormalalive %in% names(normalgenenames)
matching_genenormal <- overallnormalalive[genenormalmatches]
earlynormalgene.match <- nsclcgene.normal[, matching_genenormal]
normalmirnanames <- setNames(rep(TRUE, length(colnames(nsclcmirna.normal))), colnames(nsclcmirna.normal))
mirnanormalmatches <- overallnormalalive %in% names(normalmirnanames)
matching_mirnanormal <- overallnormalalive[mirnanormalmatches]
earlynormalmirna.match <- nsclcmirna.normal[, matching_mirnanormal]


# View the combined data frame
#print(combined_df)
#mirnansclc <- cbind(mirnaluad, mirnansclc)
#NSCLC mRNA 
genensclc <- cbind(geneluad, genensclc)


#normal gene alive 
genensclc.normal <- gsub("\\.", "-", genensclc.normal)
genenormalalive <- intersect(genensclc.normal, alive.digit)

#tumor gene alive
genesnsclc.tumor <- gsub("\\.", "-", genensclc.tumor)
genetumoralive <- intersect(genesnsclc.tumor, alive.digit)

#Obtain only alive patients
