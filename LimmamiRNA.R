# Limma #
View(rawCountresult)
dgelimma <- DGEList(counts = rawCountresult,
                    genes = rownames(rawCountresult),
                    group = sampleInfo$tumor)

# remove miRNAs that are lowly expressed
keep <- filterByExpr(dgelimma, group = dgelimma$samples$group)
dgelimma <- dgelimma[keep, , keep.lib.sizes = FALSE]

# compute scaling factors to convert observed library sizes into effective sizes
dgelimma <- calcNormFactors(dgelimma)

View(rawCountresult)
# miRNA design matrix
design <- model.matrix(~ tumor, data = sampleInfo)
rownames(design) <- colnames(dgelimma)
View(design)

# transform RNA-seq data for linear modeling
par(mar = c(1, 1, 1, 1))
dgevoom <- voom(dgelimma, design, plot = T)

# fitting the model
dgevoom.fit <- lmFit(dgevoom, design)

# perform empirical bayes smoothing standard errors
dgevoom.fit <- eBayes(dgevoom.fit)
plotSA(dgevoom.fit, main = "Final model: Mean-variance trend")

NS.resFiltlimma_miRNA <- topTable(dgevoom.fit, n = nrow(dgevoom.fit$genes),
                                  adjust.method = "BH", sort.by = "logFC", p.value = 0.01)
View(NS.resFiltlimma_miRNA)
filtlimma <- HTSFilter(NS.resFiltlimma_miRNA)$filteredData

# Downregulated
sigDownReg <- NS.resFiltlimma_miRNA[NS.resFiltlimma_miRNA$P.Value<(0.05),]
sigDownReg <- sigDownReg[sigDownReg$adj.P.Val<(0.05),]
sigDownReg <- sigDownReg[sigDownReg$logFC<(-1.5),]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
View(sigDownReg)

# Upregulated
sigUpReg <- NS.resFiltlimma_miRNA[NS.resFiltlimma_miRNA$P.Value<(0.05),]
sigUpReg <- sigUpReg[sigUpReg$adj.P.Val<(0.05),]
sigUpReg <- sigUpReg[sigUpReg$logFC>1.5,]
sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing = TRUE),]
sigUpReg <- sigUpReg[1:10, ]
View(sigUpReg)

# extract the differential expression genes from expression dataset
earlycancer_alive <- t(earlycancer_alive)
l.down.express_miRNA <- earlycancer_alive[rownames(sigDownReg),]
l.up.express_miRNA <- earlycancer_alive[rownames(sigUpReg),]
l.down.express_miRNA <- l.down.express_miRNA[rowMeans(l.down.express_miRNA == 0) < 0.4, ]
l.up.express_miRNA <- l.up.express_miRNA[rowMeans(l.up.express_miRNA == 0) < 0.4, ]
downmiRNA <- t(l.down.express_miRNA)
View(downmiRNA)
upmiRNA <- t(l.up.express_miRNA)
View(upmiRNA)
updownmiRNA <- cbind(upmiRNA, downmiRNA)
View(updownmiRNA)
write.csv(upmiRNA, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/upmiRNAlimmafix.csv")
write.csv(downmiRNA, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/downmiRNAlimmafix.csv")
write.csv(updownmiRNA, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/updownmiRNAlimmafix.csv")

#edgeR
# edgeR #

library(limma)
library(edgeR)

dgeFull <- DGEList(rawCountresult, 
                   group = sampleInfo$tumor)
dgeFull <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ], 
                   group=dgeFull$samples$group)

# Remove low expression miRNA 
keep <- filterByExpr(y = dgeFull)
dgeFull <- dgeFull[keep, , keep.lib.sizes=FALSE]
keep <-rowSums(cpm(dgeFull)>1) >= 2
dgeFull <- dgeFull[keep,]

# Estimate Normalization factor
dgeFull <- calcNormFactors(dgeFull, method = "TMM")

# Create model design matrix
design <- model.matrix(~ tumor, data = sampleInfo)
rownames(design) <- colnames(dgeFull)
View(design)

# Estimate common and tagwise dispersion
dgeFull <- estimateGLMCommonDisp(dgeFull, design=design)
dgeFull <- estimateGLMTrendedDisp(dgeFull, design=design)
dgeFull <- estimateGLMTagwiseDisp(dgeFull, design=design)
# Perform an exact test for the difference in expression between two conditions
dgeTest <- glmFit(dgeFull, design = design)
dgeTest <- glmLRT(dgeTest)
dgeTest
dgeTest <- exactTest(dgeFull)

# Extract the most differentially expressed genes
filtTMM <- HTSFilter(dgeTest)$filteredData
NS.resFilt_miRNA <- topTags(dgeTest, n=nrow(dgeTest$table),
                            sort.by = "logFC", p.value = 0.05)

# Downregulated
sigDownRegE <- NS.resFilt_miRNA$table[NS.resFilt_miRNA$table$FDR<0.05,]
sigDownRegE <- sigDownRegE[sigDownRegE$logFC<(-1.5),]
sigDownRegE <- sigDownRegE[sigDownRegE$PValue<(0.05),]
sigDownRegE <- sigDownRegE[order(sigDownRegE$logFC),]
View(sigDownRegE)

# Upregulated
sigUpRegE <- NS.resFilt_miRNA$table[NS.resFilt_miRNA$table$FDR<0.05,]
sigUpRegE <- sigUpRegE[sigUpRegE$logFC>1.5,]
sigUpRegE <- sigUpRegE[sigUpRegE$PValue<(0.05),]
sigUpRegE <- sigUpRegE[order(sigUpRegE$logFC, decreasing = TRUE),]
View(sigUpRegE)

# extract the differential expression genes from expression dataset
down.express_miRNA <- earlycancer_alive[rownames(sigDownRegE),]
up.express_miRNA <- earlycancer_alive[rownames(sigUpRegE),]
down.express_miRNA <- down.express_miRNA[rowMeans(down.express_miRNA == 0) < 0.4, ]
up.express_miRNA <- up.express_miRNA[rowMeans(up.express_miRNA == 0) < 0.4, ]
View(down.express_miRNA)
View(up.express_miRNA)
--------
write.csv(upmiRNA, )
#Merge Upregulated
rownames(upmiRNA) <- upmiRNA$Patient
View(upmiRNA)
upmiRNA <- as.data.frame(upmiRNA)
View(up.express_miRNA)
up.express_miRNA <- as.data.frame(t(up.express_miRNA))
rownames(up.express_miRNA) <- up.express_miRNA$Patient
View(up.express_miRNA)
mir.match <- intersect(rownames(upmiRNA), colnames(up.express_miRNA))
mir.match
cor.match <- as.data.frame(cor.match)
cor.match
View(upmiRNA)
mirna.match <- cor.match [, df1]
limmaedgeRmiRNA <- merge(upmiRNA, up.express_miRNA, by = "Patient")
View()
#Merge Downregulated
  
l.NS.down_miRNA <- rownames(l.down.express_miRNA)
l.NS.up_miRNA <- rownames(l.up.express_miRNA)

## overlap limma and edgeR ##
set1_NS_miRNA <- as.character(NS.down_miRNA) 
set2_NS_miRNA <- as.character(l.NS.down_miRNA)
set3_NS_miRNA <- as.character(NS.up_miRNA)
set4_NS_miRNA <- as.character(l.NS.up_miRNA)

NS.overlapping_down_miRNA <- intersect(set1_NS_miRNA, set2_NS_miRNA)
print(NS.overlapping_down_miRNA)
NS.overlapping_up_miRNA <- intersect(set3_NS_miRNA, set4_NS_miRNA)
print(NS.overlapping_up_miRNA)