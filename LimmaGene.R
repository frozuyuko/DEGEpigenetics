# Limma #
View(rawCountresult)
rawCountresult <- t(rawCountresult)
dgelimma <- DGEList(counts = rawCountresult,
                    genes = rownames(rawCountresult),
                    group = sampleinformation$Condition)

# remove miRNAs that are lowly expressed
keep <- filterByExpr(dgelimma, group = dgelimma$samples$group)
dgelimma <- dgelimma[keep, , keep.lib.sizes = FALSE]

# compute scaling factors to convert observed library sizes into effective sizes
dgelimma <- calcNormFactors(dgelimma)

View(rawCountresult)
# miRNA design matrix
design <- model.matrix(~ Condition, data = sampleinformation)
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

# Downregulated
sigDownReg <- NS.resFiltlimma_miRNA[NS.resFiltlimma_miRNA$P.Value<(0.01),]
sigDownReg <- sigDownReg[sigDownReg$adj.P.Val<(0.01),]
sigDownReg <- sigDownReg[sigDownReg$logFC<(-1.5),]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
sigDownReg <- sigDownReg[1:10, ]
View(sigDownReg)

# Upregulated
sigUpReg <- NS.resFiltlimma_miRNA[NS.resFiltlimma_miRNA$P.Value<(0.01),]
sigUpReg <- sigUpReg[sigUpReg$adj.P.Val<(0.01),]
sigUpReg <- sigUpReg[sigUpReg$logFC>1.5,]
sigUpReg <- sigUpReg[order(sigUpReg$logFC, decreasing = TRUE),]
sigUpReg <- sigUpReg[1:10, ]
View(sigUpReg)

# extract the differential expression genes from expression dataset
View(rawCountresult)
earlycanceralive <- t(rawCountresult)
geneup <- sigUpReg$genes
genedown <-sigDownReg$genes
rownames(sigUpReg) <- sigUpReg$genes
rownames(sigDownReg) <- sigDownReg$genes
l.down.express_gene <- rawCountresult[rownames(sigDownReg),]
l.up.express_gene <- rawCountresult[rownames(sigUpReg),]
l.down.express_gene <- l.down.express_miRNA[rowMeans(l.down.express_miRNA == 0) < 0.4, ]
l.up.express_gene <- l.up.express_miRNA[rowMeans(l.up.express_miRNA == 0) < 0.4, ]
downgene <- t(l.down.express_gene)
View(downgene)
upgene <- t(l.up.express_gene)
View(upgene)
sigUpDown <- cbind(upgene, downgene)
View(sigUpDown)
#upgene
write.csv(upgene, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/upgenelimma.csv")
write.csv(downgene, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/downgenelimma.csv")
write.csv(sigUpDown, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/updowngenelimma.csv")


edgeR #

library(limma)
library(edgeR)

dgeFull <- DGEList(rawCountresult, 
                   group = sampleinformation$Condition)
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
design <- model.matrix(~ Condition, data = sampleinformation)
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
NS.resFilt_gene<- topTags(dgeTest, n=nrow(dgeTest$table),
                            sort.by = "logFC", p.value = 0.05)

# Downregulated
downsig <- NS.resFilt_gene$table[NS.resFilt_gene$table$FDR<0.05,]
downsig <- downsig[downsig$logFC<(-1.5),]
downsig <- downsig[downsig$PValue<(0.05),]
downsig <- downsig[order(downsig$logFC),]
View(downsig)
# Upregulated
upsig <- NS.resFilt_gene$table[NS.resFilt_gene$table$FDR<0.05,]
upsig <- upsig[upsig$logFC>1.5,]
upsig <- upsig[upsig$PValue<(0.05),]
upsig <- upsig[order(upsig$logFC, decreasing = TRUE),]
View(upsig)

# extract the differential expression genes from expression dataset
genedata <- t(rawCountresult)
down.express_gene <- rawCountresult[rownames(downsig),]
up.express_gene <- rawCountresult[rownames(upsig),]
down.express_gene <- down.express_gene[rowMeans(down.express_gene == 0) < 0.4, ]
up.express_gene<- up.express_gene[rowMeans(up.express_gene == 0) < 0.4, ]
View(down.express_miRNA)
View(up.express_miRNA)
#Problem hereee::: 
alive_tumor <- tumor.primary[, earlytumor.match]
earlytumor.match
alive_tumor <- t(alive_tumor)
View(alive_tumor)
down.express <- alive_tumor[rownames(downgene),]
up.express <- alive_tumor[rownames(sigUpReg),]
ab.express <- alive_tumor[rownames(sigDEG),]

down.express <- down.express[apply(down.express == 0, 1, sum) <= 150,]
up.express <- up.express[apply(up.express == 0, 1, sum) <= 150,]
ab.express <- ab.express[apply(ab.express == 0, 1, sum) <= 150,]

write.csv()
