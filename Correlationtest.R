df1 <- read.csv("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/updownmiRNAedgeRlimma.csv")
df2 <- read.csv("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/updowngeneedgeRlimma.csv")
rownames(df1) <- df1$Patient
View(df1)
View(df2)
cor.match <- intersect(rownames(df1), df2$Patient)
cor.match <- as.data.frame(cor.match)
cor.match
View(df1)
mirna.match <- cor.match [, df1]
intersection <- merge(df1, df2, by = "Patient")
intersection
View(intersection)
write.csv(intersection, file = 
            "C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/cortest.csv")

#correlation test miRNA up
#miR210 to gene
mir210 <- cor(intersection$hsa.mir.210, intersection$ITLN1)
mir210neg <- mir210 * sign(mir210)
mir210neg
cor.test(intersection$hsa.mir.210, intersection$SFTPC)
cor.test(intersection$hsa.mir.210, intersection$CLDN18)
cor.test(intersection$hsa.mir.210, intersection$LGI3)
cor.test(intersection$hsa.mir.210, intersection$FABP4)
cor.test(intersection$hsa.mir.210, intersection$SERTM1)
#mir9-2
cor.test(intersection$hsa.mir.9.2, intersection$CLDN18)
cor.test(intersection$hsa.mir.9.2, intersection$SERTM1)
#mir9-1
cor.test(intersection$hsa.mir.9.1, intersection$CLDN18)
cor.test(intersection$hsa.mir.9.1, intersection$SERTM1)
#mir708
cor.test(intersection$hsa.mir.708, intersection$ANKRD1)
cor.test(intersection$hsa.mir.708, intersection$SERTM1)

cor.test(intersection$hsa.mir.196a.1, intersection$FABP4)
cor.test(intersection$hsa.mir.153.2, intersection$CLDN18)
cor.test(intersection$hsa.mir.96, intersection$CLDN18)
cor.test(intersection$hsa.mir.3607, intersection$SERTM1)

#miR135b
cor.test(intersection$hsa.mir.135b, intersection$ANKRD1)
cor.test(intersection$hsa.mir.135b, intersection$CD300LG)
cor.test(intersection$hsa.mir.135b, intersection$CLDN18)

#mirnadowngeneup
cor.test(intersection$hsa.mir.486, intersection$FAM83A)

cor.test(intersection$hsa.mir.144, intersection$FGB)
cor.test(intersection$hsa.mir.144, intersection$ABCA12)
cor.test(intersection$hsa.mir.143, intersection$MMP13)
cor.test(intersection$hsa.mir.143, intersection$ALB)
cor.test(intersection$hsa.mir.143, intersection$REG4)
cor.test(intersection$hsa.mir.143, intersection$SPINK1)
cor.test(intersection$hsa.mir.30a, intersection$ALB)
cor.test(intersection$hsa.mir.30a, intersection$ABCA12)
cor.test(intersection$hsa.mir.30a, intersection$SPINK1)
cor.test(intersection$hsa.mir.139, intersection$COL11A1)
cor.test(intersection$hsa.let.7c, intersection$FGB)
cor.test(intersection$hsa.mir.133a.1, intersection$SPINK1)
cor.test(intersection$hsa.mir.4732, intersection$REG4)
