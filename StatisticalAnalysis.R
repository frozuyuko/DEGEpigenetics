stats_updownnormal <- read.csv("C:/Users/ASUS/Downloads/TrainEpigenetics/TCGA-Assembler/top11updownregulated.csv", header = TRUE)
# Assuming you have a data frame named df
stats_updownnormal <- stats_updownnormal[-1]
View(stats_updownnormal)
stats_updownnormal <- within(stats_updownnormal, {
  colClasses <- c("numeric", "numeric", "numeric", "numeric", "numeric", 
                  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "factor")
})
colnames(stats_updownnormal)[1]
miR1 <- stats_updownnormal$hsa.mir.1269a
miR2 <- stats_updownnormal$hsa.mir.143
miR3 <- stats_updownnormal$hsa.mir.144
miR4 <- stats_updownnormal$hsa.mir.184
miR5 <- stats_updownnormal$hsa.mir.196a.1
mir6 <- stats_updownnormal$hsa.mir.196a.2
mir7 <- stats_updownnormal$hsa.mir.451a
mir8 <- stats_updownnormal$hsa.mir.4732
mir9 <- stats_updownnormal$hsa.mir.486
mir10 <- stats_updownnormal$hsa.mir.9.1
mir11 <- stats_updownnormal$hsa.mir.9.2
stats_updownnormal$status <- ifelse(stats_updownnormal$status == "cancer", 1, 
                                    ifelse(stats_updownnormal$status == "normal", 
                                           0, stats_updownnormal$status))
stats_updownnormal$status <- as.numeric(stats_updownnormal$status)
stats_updownnormal[, 1:11] <- lapply(stats_updownnormal[, 1:11], round, digits = 2)
cor.miR1R <- cor.test(miR1, stats_updownnormal$status, method = "spearman")
cor.miR2R <- cor.test(miR2, stats_updownnormal$status, method = "spearman")
cor.miR3R <- cor.test(miR3, stats_updownnormal$status, method = "spearman")
cor.miR4R <- cor.test(miR4, stats_updownnormal$status, method = "spearman")
cor.miR5R <- cor.test(miR5, stats_updownnormal$status, method = "spearman")
cor.miR6R <- cor.test(miR6, stats_updownnormal$status, method = "spearman")
cor.miR7R <- cor.test(miR7, stats_updownnormal$status, method = "spearman")
cor.miR8R <- cor.test(stats_updownnormal$hsa.mir.4732, stats_updownnormal$status, method = "spearman")
cor.miR9R <- cor.test(stats_updownnormal$hsa.mir.486, stats_updownnormal$status, method = "spearman")
cor.miR10R <- cor.test(stats_updownnormal$hsa.mir.9.1, stats_updownnormal$status, method = "spearman")
cor.miR11R <- cor.test(stats_updownnormal$hsa.mir.9.2, stats_updownnormal$status, method = "spearman")
print(cor.miR1R)
cor.miR1s <- c("spearman","spearman","spearman")
cor.miR1 <- cor(miR1, stats_updownnormal$status)
cor.test()
--- ANOVA
miR1anova<-aov(miR1~status, data=stats_updownnormal)
miR2anova<-aov(miR2~status, data=stats_updownnormal)
miR3anova<-aov(miR3~status, data=stats_updownnormal)
miR4anova<-aov(miR4~status, data=stats_updownnormal)
miR5anova<-aov(miR5~status, data=stats_updownnormal)
miR6anova<-aov(miR6~status, data=stats_updownnormal)
miR7anova<-aov(miR7~status, data=stats_updownnormal)
miR8anova<-aov(stats_updownnormal$hsa.mir.4732~status, data=stats_updownnormal)
miR9anova<-aov(stats_updownnormal$hsa.mir.486~status, data=stats_updownnormal)
miR10anova<-aov(stats_updownnormal$hsa.mir.9.1~status, data=stats_updownnormal)
miR11anova<-aov(stats_updownnormal$hsa.mir.9.2~status, data=stats_updownnormal)
anova(miR1anova)
summary(miR2anova)
summary(miR3anova)
summary(mir11anova)

