#Deseq 2 for Genes
library( "DESeq2" )
library(ggplot2)
#limma for miRNA
library(dplyr)
library(limma)
library(edgeR)
sampleInfo <- rbind(early.match, dat_normal)

deMFull <- DGEList(rawCountresult, group=sampleInfo$tumor)
deMFull <- DGEList(deMFull$counts[apply(deMFull$counts, 1, sum) != 0, ],
                   group=deMFull$samples$group)
View(deMFull)
read.delim()
datasetmirna <- as.data.frame((deMFull$counts))
View(datasetmirna)
d0 <- DGEList(datasetmirna)
d0 <- calcNormFactors(d0)
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d) # number of genes left
snames <- colnames(datasetmirna) # Sample names
cultivar <- substr(snames, 1, nchar(snames) - 2) 
time <- substr(snames, nchar(snames) - 1, nchar(snames) - 1)
cultivar
time
group <- interaction(cultivar, time)
group
mm <- model.matrix(~0 + cultivar)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
fit

resFilt <- topTags(coef(fit), n=nrow(deMFull$table), sort.by = "logFC", p.value = 0.01)





tmp
head(coef(fit))
View(coef(fit))
top.table <- topTable(coef(fit), sort.by = "P", n = Inf)
head(top.table, 20)
tmp <- eBayes(fit)




contr <- makeContrasts(hsa-let-7a-1 ~ hsa-mir-99b, levels = colnames(coef(fit)))
contr
tmp <- contrasts.fit(fit, contr)



# Extract the top genes or miRNAs based on differential expression
top_genes <- topTable(fit, coef = 1, number = n, adjust.method = "BH", p.value = 0.01)

# Extract the expression matrix for the top genes
top_expression <- expression_matrix[row.names(expression_matrix) %in% top_genes$ID, ]

# Calculate distances based on expression values
distances <- dist(t(top_expression))

# Perform MDS
mds <- cmdscale(distances, k = 2)  # Set k to the desired number of dimensions (e.g., 2)

# Plot MDS
plotMDS(mds, labels = colnames(expression_matrix), col = as.numeric(group))




plotMDS(d, col = as.numeric(snames))
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
tmp <- voom(d0, mm, plot = T)
fit <- lmFit(y, mm)
head(coef(fit))



readmirna <- read.delim(datasetmirna, row.names = 1)

redeMFull$samples
deMFull <- estimateCommonDisp(deMFull)
deMFull <- estimateTagwiseDisp(deMFull)
deMTest <- exactTest(deMFull)

deMTest <- t(deMTest)
deMTest

demframe <- as.data.frame(deMTest)
logCPM <- cpm(demframe, log=TRUE, prior.count=3)
logCPM
design <- deMFull$samples$group




fit <- lmFit(logCPM, group=deMFull$samples$group)
fit <- eBayes(fit, trend=TRUE)
fit <- treat(fit, lfc=log2(1.5), trend=TRUE)
View(fit)
topTable(fit, coef=ncol(group=deMFull$samples$group))



View(rawCountresult)
counts <- read.delim(rawCountresult, group=sampleInfo$tumor)
View(sampleInfo)
head(counts)