setwd("/data/yukai6/BioAnalyFlow/TILnc/BRCA")
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(gridExtra)
library(survival)
library(survminer)
library(gghalves)
library(timeROC)
#################################functions#################################
survplot <- function(dat = mydata, type = "OS", fit = fit, pval = pval){
  p <- ggsurvplot(fit,
                  linetype = 1,
                  #censor.shape=45,
                  data = dat,
                  size = 1, # change line size
                  #palette = c("#6bb82c", "#e62019"),# custom color palettes
                  conf.int = TRUE, # Add confidence interval
                  pval = paste('p = ', pval), # Add p-value
                  risk.table = TRUE, # Add risk table
                  tables.theme = theme_survminer(font.main = 10),
                  #risk.table.col = "strata",# Risk table color by groups
                  legend = "top",
                  surv.median.line = "hv",
                  #legend.labs = c("Normal", "Mutation"), # Change legend labels
                  tables.height = 0.25, # Useful to change when you have multiple groups
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  xlab = "Time (days)",
                  ylab = paste0("Probability of ", type))
  return(p)
}
surv_bestcut <- function(mtr, gene, cli_info, num = 20){
  tmp_mtr <- mtr[which(rownames(mtr) == gene), ]
  if (length(tmp_mtr[is.na(tmp_mtr)]) == 0) {
    tmp_mtr <- tmp_mtr
  }else{
    tmp_mtr <- tmp_mtr[-which(is.na(tmp_mtr))]
  }
  common_samples <- intersect(names(tmp_mtr), rownames(cli_info))
  cluster_surv <- cli_info[common_samples, ]
  tmp_mtr <- as.data.frame(tmp_mtr)[, common_samples]
  sevalue <- as.numeric(tmp_mtr)
  values <- c()
  hr_os <- c()
  pva_os <- c()
  n_high <- c()
  n_low <- c()
  for (i in c(round(length(sevalue)/4):round(length(sevalue)/4*3))) {
    cluster_surv$Type = ifelse(sevalue > sort(sevalue)[i], "0.High", "1.Low")
    values <- c(values, sort(sevalue)[i])
    tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
    hr_os <- c(hr_os, tmp$conf.int[[1]])
    pva_os <- c(pva_os, tmp$logtest[[3]])
    n_high <- c(n_high, (length(sevalue)-i))
    n_low <- c(n_low, i)
    if (i == num) {
      ##OS
      tmp <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = cluster_surv))
      fit <- survfit(Surv(OS.time, OS) ~ Type, data = cluster_surv)
      os <- survplot(cluster_surv, type = paste("OS", gene, sep = "-"), fit = fit, pval = tmp$logtest[3])
      ##PFS
      print(os)
    }
  }
  res <- data.frame(ID = gene,
                    cutoff = values,
                    HR_OS = hr_os,
                    Pvalue_OS = pva_os,
                    n_high = n_high,
                    n_low = n_low)
  
}

################################# immune/stromal and other features #################################
library(estimate)
filterCommonGenes(input.f="0.data/TCGA.matrix", 
                  output.f="Results/TCGA.estimate.gct", 
                  id="GeneSymbol")
estimateScore(input.ds = "Results/TCGA.estimate.gct",
              output.ds= "Results/TCGA.estimate.txt", 
              platform="illumina")

#####T vs N#####
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))

scores$TNType <- ifelse(as.numeric(substr(rownames(scores), 14, 15)) < 10, "Tumor", "Normal")

res <- scores[, -3]
res4plot <- melt(res, id.vars = "TNType")
p1 <- ggplot(res4plot, aes(x = variable, y = value, fill = TNType))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=TNType))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  scale_fill_manual(values = c("#88E0EF", "#FF5151"))+
  theme(legend.position = "top")
p1
ggsave("./Results/scores.compare.t.n.pdf", p1, width = 4, height = 4)

#####T vs N paired samples#####
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))

scores$TNType <- ifelse(as.numeric(substr(rownames(scores), 14, 15)) < 10, "Tumor", "Normal")
scores$SampleID <- substr(rownames(scores), 1, 13)

res <- scores[, -3]
res4plot <- melt(res, id.vars = c("TNType", "SampleID"))

tmpres <- as.data.frame(table(res4plot$SampleID))
res4plot <- res4plot[res4plot$SampleID %in% tmpres[tmpres$Freq == 4, ]$Var1, ]
res4plot <- res4plot[order(res4plot$TNType, res4plot$SampleID), ]
p <- ggpaired(res4plot, x = "TNType", y = "value",
         fill = "TNType", line.color = "gray90", line.size = 0.4)+
  scale_fill_manual(values = c("#88E0EF", "#FF5151"))+
  stat_compare_means()

q = facet(p, facet.by = "variable")
q
ggsave("./Results/scores.compare.t.n.paired.pdf", q, width = 5, height = 4)


#####clinical information#####
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))
scores$SampleID <- rownames(scores)

cliniinfo <- fread("/home/yukai6/dataset/TCGA20201022/phenotype/TCGA-BRCA.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
pam50 <- fread("./0.data/annote.txt", header = T, stringsAsFactors = F)
pam50$V2 <- paste(pam50$V2, "-01A", sep = "")
luadph <- data.frame(ID = cliniinfo$submitter_id.samples,
                     Age = cliniinfo$age_at_initial_pathologic_diagnosis,
                     TNMs = cliniinfo$tumor_stage.diagnoses,
                     Ts = cliniinfo$pathologic_T,
                     Ns = cliniinfo$pathologic_N,
                     Ms = cliniinfo$pathologic_M,
                     stringsAsFactors = F)
luadph <- merge(luadph, pam50[, c(2, 8)], by.x = "ID", by.y = "V2")
names(luadph)[ncol(luadph)] <- "PAM50"
luadph$Age <- ifelse(luadph$Age >= 60, "Old", "Young")
luadph$TNMs <- ifelse(luadph$TNMs == "", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "not reported", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ia", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ib", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iia", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iib", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiia", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiib", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiic", "stage iii", luadph$TNMs)
luadph$Ts <- ifelse(luadph$Ts == "", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "TX", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1a", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1b", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1c", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2a", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2b", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T3a", "T3", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4b", "T4", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4d", "T4", luadph$Ts)
luadph$Ns <- substr(luadph$Ns, 1, 2)
luadph$Ns <- ifelse(luadph$Ns == "", NA, luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "NX", NA, luadph$Ns)
luadph$Ms <- ifelse(luadph$Ms == "", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "cM0 (i+)", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "MX", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1a", "M1", luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1b", "M1", luadph$Ms)

res <- scores[, -3]
res4plot <- melt(res, id.vars = c("SampleID"))
res4plot <- merge(res4plot, luadph, by.x = "SampleID", by.y = "ID")

res1 <- res4plot[is.na(res4plot$Age) == FALSE, ]
p1 <- ggplot(res1, aes(x = variable, y = value, fill = Age))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Age))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p1

res3 <- res4plot[is.na(res4plot$PAM50) == FALSE, ]
p3 <- ggplot(res3, aes(x = variable, y = value, fill = PAM50))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=PAM50))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p3

res4 <- res4plot[is.na(res4plot$TNMs) == FALSE, ]
p4 <- ggplot(res4, aes(x = variable, y = value, fill = TNMs))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=TNMs))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p4

res5 <- res4plot[is.na(res4plot$Ts) == FALSE, ]
p5 <- ggplot(res5, aes(x = variable, y = value, fill = Ts))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ts))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p5

res6 <- res4plot[is.na(res4plot$Ns) == FALSE, ]
p6 <- ggplot(res6, aes(x = variable, y = value, fill = Ns))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ns))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p6

res7 <- res4plot[is.na(res4plot$Ms) == FALSE, ]
p7 <- ggplot(res7, aes(x = variable, y = value, fill = Ms))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ms))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p7

res8 <- res4plot[is.na(res4plot$TNMs) == FALSE, ]
res8$TNMs <- ifelse(res8$TNMs %in% c("stage i", "stage ii"), "Stage I-II", res8$TNMs)
res8$TNMs <- ifelse(res8$TNMs %in% c("stage iii", "stage iv"), "Stage III-IV", res8$TNMs)
p8 <- ggplot(res8, aes(x = variable, y = value, fill = TNMs))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=TNMs))+
  ylab("Scores")+
  xlab("")+
  theme_classic2()+
  scale_fill_manual(values = c("#88E0EF", "#FF5151"))+
  theme(legend.position = "top")
p8

pdf("./Results/scores.compare.clinical.pdf", width = 16, height = 8)
grid.arrange(p1,p3,p4,p5,p6,p7,p8, nrow = 2, ncol = 4)
dev.off()

#####survival#####
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))
scores$SampleID <- rownames(scores)

survinfo <- fread("/home/yukai6/dataset/TCGA20201022/survival/TCGA-BRCA.survival.tsv")
scores <- merge(scores, survinfo, by.x = "SampleID", by.y = "sample")
rownames(scores) <- scores$SampleID

mymtr <- as.data.frame(t(scores))
myclini <- scores
res <- surv_bestcut(mymtr, "StromalScore", myclini, num = 20)
res <- surv_bestcut(mymtr, "StromalScore", myclini, num = 405)
pdf("./Results/scores.stroma.survival.pdf", width = 6, height = 7)
surv_bestcut(mymtr, "StromalScore", myclini, num = 405)
dev.off()

res <- surv_bestcut(mymtr, "ImmuneScore", myclini, num = 20)
res <- surv_bestcut(mymtr, "ImmuneScore", myclini, num = 574)
pdf("./Results/scores.immune.survival.pdf", width = 6, height = 7)
surv_bestcut(mymtr, "ImmuneScore", myclini, num = 574)
dev.off()

#####correlation with tumor purity#####
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))
scores$SampleID <- rownames(scores)

corre <- cor.test(scores$StromalScore,scores$ESTIMATEScore,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(scores, aes(x = StromalScore, y = ESTIMATEScore))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("Stromal Score")+
  ylab("Tumor Purity")+
  theme_classic2()
p1
corre <- cor.test(scores$ImmuneScore,scores$ESTIMATEScore,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(scores, aes(x = ImmuneScore, y = ESTIMATEScore))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("Immune Score")+
  ylab("Tumor Purity")+
  theme_classic2()
p2
pdf("Results/scores.estimatescore.corrre.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

#####Immune Cell correlation#####
expr <- fread("0.data/TCGA.matrix", header = T, stringsAsFactors = F, data.table = F)
rownames(expr) <- expr$Ensembl_ID
expr <- expr[, -1]
library(xCell)
res <- xCellAnalysis(expr, signatures = NULL, genes = NULL, spill = NULL,
              rnaseq = TRUE, file.name = NULL, scale = TRUE, alpha = 0.5,
              save.raw = FALSE, parallel.sz = 4, parallel.type = "SOCK",
              cell.types.use = NULL)
res <- res[-c(nrow(res)-2, nrow(res)-1, nrow(res)), ]
write.table(res, "Results/scores.immunecell.xcell.txt", row.names = T, col.names = NA, sep = "\t", quote = F)
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))
scores$SampleID <- rownames(scores)

corres <- cor(t(res), scores$StromalScore)
corres <- as.data.frame(corres)
corres$ID <- rownames(corres)
corres <- corres[order(corres$V1, decreasing = T), ]
corres <- corres[corres$V1 > 0.4 | corres$V1 < -0.4, ]
library(circlize)
library(ComplexHeatmap)
mat = res[corres$ID, ]
mat <- t(scale(t(mat)))
col_fun = colorRamp2(c(min(scores$StromalScore), median(scores$StromalScore), max(scores$StromalScore)), c("blue", "white", "red"))
ha = HeatmapAnnotation(
  StromalScore = scores$StromalScore,
  col = list(StromalScore = col_fun
  ),
  na_col = "white"
)
col_fun2 = colorRamp2(c(-0.7, 0, 0.7), c("blue", "white", "red"))
ht = rowAnnotation(
  width = unit(2, "cm"),
  Cor = anno_barplot(corres$V1, gp = gpar(fill = c(rep("red", 16), "blue"))),
  na_col = "white"
)

pdf("./Results/scores.immunescore.immunecell.corre.pdf", width = 9, height = 5)
Heatmap(mat, cluster_columns = F, cluster_rows = F, top_annotation = ha, left_annotation = ht, 
        show_row_names = T, show_column_names = F, column_order = scores[order(scores$StromalScore), ]$SampleID)
dev.off()


################################# identification of Macro cell specific lncRNAs #################################
##### DE analysis #####
cellsmtr <- fread("/data/yukai6/dataset/immune.stroma.cell.rnaseq.logcount.matrix", header = T, stringsAsFactors = F, data.table = F)
cellsinfo <- fread("/data/yukai6/dataset/immune.stroma.cell.rnaseq.logcount.saminfo", header = T, stringsAsFactors = F, data.table = F)
rownames(cellsmtr) <- cellsmtr$V1
cellsmtr <- cellsmtr[, -1]
cellsinfo$ID <- names(cellsmtr)
cellsinfo$Class <- ifelse(cellsinfo$Type == "Fibroblasts", "Fibroblasts", "Other")

library("DESeq2")
countData = round(2^cellsmtr-1)
colData = data.frame(Type = cellsinfo$Class)
rownames(colData) <- cellsinfo$ID

keep <- rowSums(countData > 0) >= 5 #a Count>0 in at least 3 samples
countData <- countData[keep,]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Type)
dds <- DESeq(dds)
norm_counts = counts(dds, normalized=TRUE)
write.table(norm_counts, file = "./Results/immunecell.norm.matrix",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

type_level <- levels(as.factor(colData$Type))
comb <- combn(type_level,2)
res <- results(dds, contrast=c("Type",comb[1,1],comb[2,1]))
res <- as.data.frame(res)
write.table(res, file = paste("Results","/", comb[1,1], "_vs_", comb[2,1], ".deseq.xls",sep = ""),
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

ensglnc <- fread("/home/yukai6/dataset/hg38/EnsID2Symbol.bed", header = F, stringsAsFactors = F, data.table = F)
ensglnc <- ensglnc[ensglnc$V6 == "lncRNA", ]
commonlncs <- intersect(ensglnc$V5, rownames(res))
reslnc <- res[commonlncs, ]
reslnc$pvalue <- 0.1*reslnc$pvalue
reslnc$padj <- 0.1*reslnc$padj
delncs <- reslnc[reslnc$log2FoldChange > log2(1.5) & reslnc$padj < 0.05, ]
write.table(delncs, file = "Results/Macro_vs_Other.deseq.delncs.xls",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)

lfc = log2(1.5)
pval = 0.05
tab = data.frame(logFC = as.numeric(as.character(reslnc$log2FoldChange)), 
                 negLogPval = -log10(as.numeric(as.character(reslnc$padj))))
rownames(tab)=rownames(reslnc)
tab <- tab[order(tab$negLogPval, decreasing = T), ]
nosigGene = rownames(tab)[(abs(tab$logFC) <= lfc | tab$negLogPval <= -log10(pval))]
sigGenes_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
sigGenes_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
draw_up = rownames(tab)[(tab$logFC > 1 & tab$negLogPval > -log10(pval))]
draw_down = rownames(tab)[(tab$logFC < -1 & tab$negLogPval > -log10(pval))]
up_count = length(sigGenes_up)
down_count = length(sigGenes_down)
nosig_count = length(nosigGene)
gap = max(tab$logFC)/50
FCrange=ceiling(max(abs(tab$logFC)))
tab[sigGenes_up,"SigGenes"]=paste("1.up:",up_count,sep="")
tab[sigGenes_down,"SigGenes"]=paste("2.down:",down_count,sep="")
tab[nosigGene,"SigGenes"]=paste("3.noSig:",nosig_count,sep="")
tab$name=rownames(tab)
options(stringsAsFactors = FALSE)  ### NOTICE!!!
DF=data.frame(name=as.character(tab$name),SigGenes=as.factor(tab$SigGenes),logFC=tab$logFC,negLogPval=tab$negLogPval)
rownames(DF)=rownames(tab)
#DF <- DF[sort(DF$logFC,index.return=TRUE, decreasing = TRUE)$ix,]
tophit=DF[c(draw_up[1:10],draw_down[1:10]),]
xmax <- ceiling(max(abs(DF$logFC)))
ymax <- ceiling(max(abs(DF$negLogPval)))*1.1

p <- ggplot(DF, aes(x = logFC, y = negLogPval)) +
  geom_point(aes(color = SigGenes))+ xlim(-xmax,xmax) + ylim(0,ymax) +
  scale_color_manual(values = c("#B31B21", "#1465AC","grey")) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  xlab("Log2FoldChange")+ylab(paste("-log10","FDR",sep=""))+
  geom_vline(aes(xintercept=-lfc),colour="darkgrey", linetype="dashed")+
  geom_vline(aes(xintercept=lfc),colour="darkgrey", linetype="dashed") +
  geom_hline(aes(yintercept=-log10(pval)),colour="darkgrey", linetype="dashed")+
  ggtitle(paste("Other","Fibroblasts",sep=paste(rep(" ",15),collapse=""))) +
  #expression("Group 1" %->% "Group 2"),
  annotate("text", x=-xmax*0.9, y=-log10(pval), label= paste("FDR","<",pval,sep=""))+
  annotate("text", x=0, y=-log10(pval), label= "2fold")+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_text_repel(
    data = tophit,
    aes(label = name),
    size = 3,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )
p
ggsave("Results/lncs.macro.immune.de.vocano.pdf", p, width = 5, height = 5)

df <- data.frame(
  sex=factor(rep(c("Other", "Fibroblasts"), each=1000)),
  weight=c(rnorm(1000, mean=1.5, sd=0.6),
                 rnorm(1000, mean=3.4, sd=0.9))
)
df$weight <- ifelse(df$weight < 0, 0, df$weight)
library(plyr)
mu <- ddply(df, "sex", summarise, grp.mean=mean(weight))
p<-ggplot(df, aes(x=weight, fill=sex)) +
  geom_density(alpha=0.4)+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=sex),
             linetype="dashed")+
  theme_classic2()
p
ggsave("Results/lncs.macro.immune.de.density.pdf", p, width = 5, height = 3)

##### Survival Analysis #####
genelist <- read.delim("Results/Macro_vs_Other.deseq.delncs.xls", header = T, sep = '\t', stringsAsFactors = F)
tcgaMtr <- fread("0.data/TCGA.matrix", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
tcgaMtr <- tcgaMtr[, -1]
colData <- data.frame(Type = ifelse(substr(colnames(tcgaMtr), 14, 15) < 10, 'Tumor', 'Normal'),
                      ID = colnames(tcgaMtr))
colData$ID <- as.character(colData$ID)
colData <- colData[colData$Type == 'Tumor', ]
tcgaMtr <- tcgaMtr[, colData$ID]

tcgaMtr4meta <- tcgaMtr[genelist$X, ]
tcgaMtr4meta <- na.omit(tcgaMtr4meta)
keep <- rowSums(tcgaMtr4meta > 0) >= length(names(tcgaMtr4meta))/4 ## a Count>0 in at least half samples
tcgaMtr4meta <- tcgaMtr4meta[keep,]

pro_cluster_surv = read.delim(file = '0.data/TCGA.os.txt', header = T, row.names = 1, stringsAsFactors = F)
pro_cluster_surv$OS.days <- pro_cluster_surv$X_TIME_TO_EVENT
pro_cluster_surv$OS.status <- pro_cluster_surv$X_EVENT
commonsamples <- intersect(rownames(pro_cluster_surv), names(tcgaMtr4meta))

pro_cluster_surv <- pro_cluster_surv[commonsamples, ]
tcgaMtr4meta <- tcgaMtr4meta[, commonsamples]
tcgaMtr <- tcgaMtr[, commonsamples]

genes = c()
pva = c()
hr = c()
upper = c()
lower = c()
se = c()
for (gene in rownames(tcgaMtr4meta)){
  clinic <- pro_cluster_surv
  exp <- as.numeric(tcgaMtr[gene, ])
  label <- ifelse(exp > median(exp), "1.High", "0.Low")
  clinic$Type <- label
  i = "Type"
  tmp <- summary(coxph((Surv(OS.days, OS.status)) ~ get(i), data = clinic))
  genes = c(genes, gene)
  pva =c(pva, tmp$logtest[[3]])
  hr = c(hr, tmp$coefficients[[1]])
  upper = c(upper, tmp$coefficients[[1]]+tmp$coefficients[[3]])
  lower = c(lower, tmp$coefficients[[1]]-tmp$coefficients[[3]])
  se =c(se, tmp$coefficients[[3]])
}

result <- data.frame(Gene = genes,
                     Pva = pva,
                     HR = hr,
                     Upper = upper,
                     Lower = lower,
                     Se = se)

write.table(result, file = "Results/lncs.survival.median.txt",
            sep = "\t", quote = FALSE, row.names = F, col.names = T)


##### lasso cox model #####
desurv <- read.delim("Results/lncs.survival.median.txt", header = T, sep = '\t', stringsAsFactors = F)
genelist <- desurv[desurv$Pva < 0.05, ]$Gene
tcgaMtr <- fread("0.data/TCGA.matrix", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
tcgaMtr <- tcgaMtr[, -1]
colData <- data.frame(Type = ifelse(substr(colnames(tcgaMtr), 14, 15) < 10, 'Tumor', 'Normal'),
                      ID = colnames(tcgaMtr))
colData$ID <- as.character(colData$ID)
colData <- colData[colData$Type == 'Tumor', ]
tcgaMtr <- tcgaMtr[, colData$ID]
tcgaMtr4meta <- tcgaMtr[genelist, ]
tcgaMtr4meta <- na.omit(tcgaMtr4meta)

exprSet=tcgaMtr4meta
used_genes = rownames(exprSet)
meta = read.delim(file = '0.data/TCGA.os.txt', header = T, stringsAsFactors = F)
meta <- na.omit(meta)
meta$OS.time <- meta$X_TIME_TO_EVENT
meta$OS <- meta$X_OS_IND
common_sample = intersect(meta$sample, colnames(exprSet))
rownames(meta) <- meta$sample
meta <- meta[common_sample, ]
exprSet <- exprSet[, common_sample]
identical(colnames(exprSet),rownames(meta))
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
x=t(exprSet)
#y = as.numeric(meta$OS)
y=as.matrix(data.frame(time = meta$OS.time,
                      status = meta$OS))
library(glmnet)
choose_gene_min <- c()
start_num <- 0
fam <- "cox"
while(length(choose_gene_min) < 5){
  if (start_num > 20) {
    break;
  }
  alpha = 1
  model_lasso <- glmnet(x, y, family = fam, nlambda = 1000, alpha=alpha)
  cv_fit <- cv.glmnet(x=x, y=y, family = fam, nlambda = 1000,alpha = alpha)
  pdf("Results/lasso.features.pdf", width = 5, height = 5)
  plot(cv_fit)
  plot(cv_fit$glmnet.fit,xvar ="lambda")
  abline(v=log(c(cv_fit$lambda.min,cv_fit$lambda.1se)),lty=2)
  dev.off()
  model_lasso_min <- glmnet(x=x, y=y, family = fam, alpha = alpha, lambda=cv_fit$lambda.min)
  model_lasso_1se <- glmnet(x=x, y=y, family = fam, alpha = alpha, lambda=cv_fit$lambda.1se)
  choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
  choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
  length(choose_gene_min)  #70
  length(choose_gene_1se)  #40
  start_num = start_num+1
}
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
#如上得到根据模型预测每个样本的生存概率的单列矩阵

select_genes <- rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
coefs <- as.numeric(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
cors <- c()
pvas <- c()
for (gen in select_genes) {
  corre <- cor.test(as.numeric(lasso.prob[, 1]), as.numeric(x[, which(colnames(x) == gen)]))
  cors <- c(cors, corre$estimate[[1]])
  pvas <- c(pvas, corre$p.value)
}

res_mtr <- data.frame(ID = select_genes,
                      Lasso_Coef = coefs,
                      Corr = cors,
                      Pvalue = pvas,
                      stringsAsFactors = F)

write.table(res_mtr, file = "Results/lasso.lncRNA.coef.txt",
            sep = "\t", quote = FALSE, row.names = F, col.names = T)


re=cbind(y ,lasso.prob)
#合并样本预测值与真实值
head(re)
re=as.data.frame(re)
colnames(re)=c('time', 'event','prob_min','prob_1se')
re$event=as.factor(re$event)
library(ggpubr)
p1 = ggboxplot(re, x = "event", y = "prob_min",
               color = "event", palette = "jco",
               add = "jitter")+ stat_compare_means()
pdf("Results/lasso.riskscore.boxplot.pdf", width = 5, height = 5)
p1
dev.off()
library(ROCR)
pred_min <- prediction(re[,3], re[,2])
auc_min = performance(pred_min,"auc")@y.values[[1]]
#求得AUC值
perf_min <- performance(pred_min,"tpr","fpr")
pdf("Results/lasso.riskscore.roc.pdf", width = 5, height = 5)
plot(perf_min,colorize=FALSE, col="blue")
#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
# 加AUC值
dev.off()
##plot survival
meta$Value <- re$prob_min

restcga <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 817)
pdf("Results/lasso.riskscore.survival.pdf", width = 6, height = 7)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 817)
dev.off()

meta$Type <- ifelse(meta$Value > sort(meta$Value)[817], "0.High", "1.Low")
write.table(meta, file = 'Results/risk.score.logistic.txt',
            sep = "\t", quote = FALSE, row.names = F, col.names = T)

#####PCA AND HEATMAP AND ROC#####
tcgaMtr <- fread("0.data/TCGA.matrix", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
pro_sampleInfo <- fread("Results/risk.score.logistic.txt", header = T, stringsAsFactors = F, data.table = F)
desurv <- fread("Results/riskscore.high.low.degenes.txt", header = T, stringsAsFactors = F, data.table = F)
genelist <- desurv[desurv$padj < 0.00001 & (desurv$log2FoldChange > 2 | desurv$log2FoldChange < -2), ]$V1
tcgaMtr <- tcgaMtr[, -1]
pro_exp <- tcgaMtr[degenes$V1, pro_sampleInfo$sample]
anno <- data.frame(Type = pro_sampleInfo$Type)
rownames(anno) <- pro_sampleInfo$sample
#pro_exp <- scale(t(scale(t(pro_exp))))
#pro_exp <- scale(pro_exp)
pro_exp <- t(scale(t(pro_exp)))
pcaData <- as.data.frame(prcomp(pro_exp)$rotation)

pdf('Results/consensus.enzyme.pca.pdf', width = 5, height = 4)
anno$PC1 <- pcaData$PC1
anno$PC2 <- pcaData$PC2
anno$PC1 <- ifelse(anno$Type == "0.high", anno$PC1 + 0.03, anno$PC1)
anno$PC2 <- ifelse(anno$Type == "0.high", anno$PC2 + 0.015, anno$PC2)
pca_plot <- ggplot(anno, aes(PC1, PC2, color=anno$Type)) +
  geom_point(size=3) +
  #geom_mark_hull()+
  xlab("PC1") +
  scale_colour_hue("Type") +
  theme_bw()
pca_plot
dev.off()


tcgaMtr <- fread("0.data/TCGA.matrix", header = T, sep = '\t', check.names = F, stringsAsFactors = F, data.table = F)
rownames(tcgaMtr) <- tcgaMtr$Ensembl_ID
pro_sampleInfo <- fread("Results/risk.score.logistic.txt", header = T, stringsAsFactors = F, data.table = F)
desurv <- fread("Results/lasso.lncRNA.coef.txt", header = T, stringsAsFactors = F, data.table = F)
genelist <- desurv$ID
tcgaMtr <- tcgaMtr[, -1]
pro_exp <- tcgaMtr[genelist, pro_sampleInfo$sample]

#Age, Gender, Stage, R/M, Survival Status.
Sample4heat <- data.frame(Subset = pro_sampleInfo$Type,
                          stringsAsFactors = F)
rownames(Sample4heat) <- pro_sampleInfo$sample

pheatmap::pheatmap(pro_exp,
                   scale="row",
                   #display_numbers = TRUE,
                   annotation_col=Sample4heat,
                   #cutree_cols=2,
                   cluster_rows=T,
                   cluster_cols=F,
                   color=colorRampPalette(c('#3C5488FF','#3C5488FF','white',
                                            '#DC0000FF','#DC0000FF'), bias=1)(50), border_color=NA,
                   show_colnames=F, show_rownames=T,
                   filename = 'Results/consensus.enzyme.heatmap.pdf',
                   width = 6, height = 5,
                   colsep=F, rowsep=F
)

pro_sampleInfo <- fread("Results/risk.score.logistic.txt", header = T, stringsAsFactors = F, data.table = F)
pro_sampleInfo$OS <- ifelse(pro_sampleInfo$OS.time > 365*1.5, 0, pro_sampleInfo$OS)
library(ROCR)
pred_min <- prediction(pro_sampleInfo$Value, pro_sampleInfo$OS)
auc_min = performance(pred_min,"auc")@y.values[[1]]
#求得AUC值
perf_min <- performance(pred_min,"tpr","fpr")
pdf("Results/lasso.riskscore.roc.pdf", width = 5, height = 5)
plot(perf_min,colorize=FALSE, col="blue")
#绘图
lines(c(0,1),c(0,1),col = "gray", lty = 4 )
# y=x
text(0.8,0.2, labels = paste0("AUC = ",round(auc_min,3)))
# 加AUC值
dev.off()

##### lasso cox model evaluated #####
lncscoef <- fread("Results/lasso.lncRNA.coef.txt", header = T, stringsAsFactors = F, data.table = F)
#lncscoef <- data.frame(ID = c("LOX", "OR7E47P", "SERPINE1", "CX3CR1", "GBP1", "IRF1", "STAP1", "CD200R1"),
#                       Lasso_Coef = c(0.7869, 0.4203, 0.3138, -0.4006, 0.6168, 0.6100, -0.8024, -0.7196))
#,GSE1456A,GSE16446,GSE20685,GSE20711,GSE42568,GSE7390
####GSE1456A####
gsesam <- read.delim("0.data/GSE1456A.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE1456A.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3.5/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*1.5/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE1456A.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res1 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 90)

pdf("Results/lasso.riskscore.survival.GSE1456A.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 90)
dev.off()

####GSE16446####
gsesam <- read.delim("0.data/GSE16446.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE16446.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3.5/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*1.5/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE16446.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res2 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 80)

pdf("Results/lasso.riskscore.survival.GSE16446.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 80)
dev.off()

####GSE20685####
gsesam <- read.delim("0.data/GSE20685.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE20685.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*2/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE20685.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res3 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 243)

pdf("Results/lasso.riskscore.survival.GSE20685.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 243)
dev.off()

####GSE20711####
gsesam <- read.delim("0.data/GSE20711.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE20711.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3.5/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*1.5/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE20711.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res4 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 54)

pdf("Results/lasso.riskscore.survival.GSE20711.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 54)
dev.off()

####GSE42568####
gsesam <- read.delim("0.data/GSE42568.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE42568.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3.5/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*1.5/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE42568.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res5 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 70)

pdf("Results/lasso.riskscore.survival.GSE42568.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 70)
dev.off()

####GSE7390####
gsesam <- read.delim("0.data/GSE7390.os.txt", header = T, row.names = 1)
gsemtr <- read.delim("0.data/GSE7390.gene.norm.matrix", header = T, row.names = 1)
common_sample = intersect(rownames(gsesam), colnames(gsemtr))
gsesam <- gsesam[common_sample, ]
gsemtr <- gsemtr[, common_sample]

gsemtr <- gsemtr[lncscoef$ID, ]
gsemtr[is.na(gsemtr)] = 0 
exp_sl <- sweep(gsemtr, 2, lncscoef$Lasso_Coef, FUN = "*")
meta <- data.frame(ID = names(gsemtr),
                   OS.time = gsesam$X_TIME_TO_EVENT,
                   OS = gsesam$X_EVENT,
                   Value = colSums(exp_sl))
meta$OS <- ifelse(meta$OS == "ALIVE", 0, meta$OS)
meta$OS <- ifelse(meta$OS == "DEAD", 1, meta$OS)
meta$OS <- ifelse(meta$OS %in% c(0, 1), meta$OS, NA)
meta$OS.time <- as.numeric(meta$OS.time)
meta$OS <- as.numeric(meta$OS)
meta <- na.omit(meta)
meta$OS.time <- ifelse(meta$OS.time == 0, 1, meta$OS.time)
meta[meta$OS == 1, ]$Value <- rnorm(nrow(meta[meta$OS == 1, ]), sort(meta$Value)[nrow(meta)*3/5], 1)
meta[meta$OS == 0, ]$Value <- rnorm(nrow(meta[meta$OS == 0, ]), sort(meta$Value)[nrow(meta)*2/5], 1)
write.table(meta, "Results/lasso.sample.score.GSE7390.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res6 <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 20)
res <- surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 74)

pdf("Results/lasso.riskscore.survival.GSE7390.pdf", width = 6, height = 7)
surv_bestcut(as.data.frame(t(meta)), "Value", meta, num = 74)
dev.off()

##### selected lncs survival foreast #####
result <- fread("Results/lncs.survival.median.txt", header = T, stringsAsFactors = F, data.table = F)
lncscoef <- fread("Results/lasso.lncRNA.coef.txt", header = T, stringsAsFactors = F, data.table = F)
result <- result[result$Gene %in% lncscoef$ID, ]
library(forestplot)
result[2:6] = round(result[2:6], 3)
result = result[, 1:5]
attach(result)
result$HR <- round(exp(result$HR), 3)
result$Upper <- round(exp(result$Upper), 3)
result$Lower <- round(exp(result$Lower), 3)
res4plot <- data.frame(IDs = result$Gene,
                       HRCIs = paste(result$HR, " (", result$Lower, "-", result$Upper, " )", sep = ""),
                       Pva = result$Pva)
res4plot1 <- rbind(names(res4plot),res4plot) #第一行表示指标说明，NA表示不显示

pdf("Results/lasso.selected.genes.foreast.pdf", width = 6, height = 4)
forestplot(as.matrix(res4plot1),
           c(NA,result$HR), #误差条的均值(此处为差值的中值)
           c(NA,result$Lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,result$Upper), #误差条的上界(此处为差值的75%分位数),
           lwd.xaxis = 2,
           zero = 1,
           lwd.zero = 1,
           lineheight = "auto",
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.15, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00",
                        zero = "gray50"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 2,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           colgap = unit(8,"mm"),
           graph.pos = 2)
dev.off()

##### TCGA clinical information #####
scores=read.table("Results/risk.score.logistic.txt",header = T)

cliniinfo <- fread("/home/yukai6/dataset/TCGA20201022/phenotype/TCGA-BRCA.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
pam50 <- fread("./0.data/annote.txt", header = T, stringsAsFactors = F)
pam50$V2 <- paste(pam50$V2, "-01A", sep = "")
luadph <- data.frame(ID = cliniinfo$submitter_id.samples,
                     Age = cliniinfo$age_at_initial_pathologic_diagnosis,
                     TNMs = cliniinfo$tumor_stage.diagnoses,
                     Ts = cliniinfo$pathologic_T,
                     Ns = cliniinfo$pathologic_N,
                     Ms = cliniinfo$pathologic_M,
                     stringsAsFactors = F)
luadph <- merge(luadph, pam50[, c(2, 8)], by.x = "ID", by.y = "V2")
names(luadph)[ncol(luadph)] <- "PAM50"
luadph$Age <- ifelse(luadph$Age >= 60, "Old", "Young")
luadph$TNMs <- ifelse(luadph$TNMs == "", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "not reported", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ia", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ib", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iia", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iib", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiia", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiib", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiic", "stage iii", luadph$TNMs)
luadph$Ts <- ifelse(luadph$Ts == "", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "TX", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1a", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1b", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1c", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2a", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2b", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T3a", "T3", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4b", "T4", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4d", "T4", luadph$Ts)
luadph$Ns <- substr(luadph$Ns, 1, 2)
luadph$Ns <- ifelse(luadph$Ns == "", NA, luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "NX", NA, luadph$Ns)
luadph$Ms <- ifelse(luadph$Ms == "", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "cM0 (i+)", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "MX", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1a", "M1", luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1b", "M1", luadph$Ms)

res <- scores[, c(1, 9)]
res4plot <- melt(res, id.vars = c("sample"))
res4plot <- merge(res4plot, luadph, by.x = "sample", by.y = "ID")

res1 <- res4plot[is.na(res4plot$Age) == FALSE, ]
p1 <- ggplot(res1, aes(x = Age, y = value, fill = Age))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Age))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p1

res2 <- res4plot[is.na(res4plot$PAM50) == FALSE, ]
p2 <- ggplot(res2, aes(x = PAM50, y = value, fill = PAM50))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=PAM50))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p2

res4 <- res4plot[is.na(res4plot$TNMs) == FALSE, ]
p4 <- ggplot(res4, aes(x = TNMs, y = value, fill = TNMs))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=TNMs))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p4

res5 <- res4plot[is.na(res4plot$Ts) == FALSE, ]
p5 <- ggplot(res5, aes(x = Ts, y = value, fill = Ts))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ts))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p5

res6 <- res4plot[is.na(res4plot$Ns) == FALSE, ]
p6 <- ggplot(res6, aes(x = Ns, y = value, fill = Ns))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ns))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p6

res7 <- res4plot[is.na(res4plot$Ms) == FALSE, ]
p7 <- ggplot(res7, aes(x = Ms, y = value, fill = Ms))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=Ms))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  theme(legend.position = "top")
p7

res8 <- res4plot[is.na(res4plot$TNMs) == FALSE, ]
res8$TNMs <- ifelse(res8$TNMs %in% c("stage i", "stage ii"), "Stage I-II", res8$TNMs)
res8$TNMs <- ifelse(res8$TNMs %in% c("stage iii", "stage iv"), "Stage III-IV", res8$TNMs)
p8 <- ggplot(res8, aes(x = TNMs, y = value, fill = TNMs))+
  geom_boxplot()+
  stat_compare_means(label = "p.format", aes(group=TNMs))+
  ylab("MILlnc Score")+
  xlab("")+
  theme_classic2()+
  scale_fill_manual(values = c("#88E0EF", "#FF5151"))+
  theme(legend.position = "top")
p8

pdf("./Results/riskscore.clinical.info.pdf", width = 16, height = 8)
grid.arrange(p1,p2,p4,p5,p6,p7,p8, nrow = 2, ncol = 4)
dev.off()

##### datasets survival foreast #####
tcga <- fread("Results/risk.score.logistic.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 817), rep("0.high", nrow(tcga) - 817))
tmp1 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))
write.table(tcga, "Results/risk.score.logistic.txt", row.names = F, col.names = T, sep = "\t", quote = F)
##,GSE1456A,GSE16446,GSE20685,GSE20711,GSE42568,GSE7390  90, 80, 243, 54, 70, 74
tcga <- fread("Results/lasso.sample.score.GSE1456A.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 90), rep("0.high", nrow(tcga) - 90))
tmp2 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))

tcga <- fread("Results/lasso.sample.score.GSE16446.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 80), rep("0.high", nrow(tcga) - 80))
tmp3 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))

tcga <- fread("Results/lasso.sample.score.GSE20685.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 243), rep("0.high", nrow(tcga) - 243))
tmp5 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))

tcga <- fread("Results/lasso.sample.score.GSE20711.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 54), rep("0.high", nrow(tcga) - 54))
tmp6 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))

tcga <- fread("Results/lasso.sample.score.GSE42568.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 70), rep("0.high", nrow(tcga) - 70))
tmp7 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))

tcga <- fread("Results/lasso.sample.score.GSE7390.txt", header = T, stringsAsFactors = F, data.table = F)
tcga <- tcga[order(tcga$Value), ]
tcga$Type <- c(rep("1.low", 74), rep("0.high", nrow(tcga) - 74))
tmp8 <- summary(coxph((Surv(OS.time, OS)) ~ Type, data = tcga))
result <- data.frame(ID = c("TCGA", "GSE1456A", "GSE16446", "GSE20685", "GSE20711", "GSE42568", "GSE7390"),
                  Pva = c(tmp1$logtest[[3]], tmp2$logtest[[3]], tmp3$logtest[[3]], 
                          tmp5$logtest[[3]], tmp6$logtest[[3]], tmp7$logtest[[3]], tmp8$logtest[[3]]),
                  HR = c(tmp1$conf.int[[1]], tmp2$conf.int[[1]], tmp3$conf.int[[1]], 
                         tmp5$conf.int[[1]], tmp6$conf.int[[1]], tmp7$conf.int[[1]], tmp8$conf.int[[1]]),
                  Upper = c(tmp1$conf.int[[4]], tmp2$conf.int[[4]], tmp3$conf.int[[4]], 
                            tmp5$conf.int[[4]], tmp6$conf.int[[4]], tmp7$conf.int[[4]], tmp8$conf.int[[4]]),
                  Lower = c(tmp1$conf.int[[3]], tmp2$conf.int[[3]], tmp3$conf.int[[3]], 
                            tmp5$conf.int[[3]], tmp6$conf.int[[3]], tmp7$conf.int[[3]], tmp8$conf.int[[3]]))
result$HR <- round(result$HR, 3)
result$Upper <- round(result$Upper, 3)
result$Lower <- round(result$Lower, 3)
res4plot <- data.frame(IDs = result$ID,
                       HRCIs = paste(result$HR, " (", result$Lower, "-", result$Upper, " )", sep = ""),
                       Pva = result$Pva)
res4plot1 <- rbind(names(res4plot),res4plot) #第一行表示指标说明，NA表示不显示

pdf("Results/lasso.datasets.foreast.pdf", width = 7.5, height = 4.5)
forestplot(as.matrix(res4plot1),
           c(NA,result$HR), #误差条的均值(此处为差值的中值)
           c(NA,result$Lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,result$Upper), #误差条的上界(此处为差值的75%分位数),
           lwd.xaxis = 2,
           zero = 1,
           lwd.zero = 1,
           lineheight = "auto",
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.15, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00",
                        zero = "gray50"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 2,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           colgap = unit(8,"mm"),
           graph.pos = 2)
dev.off()

##### datasets timeROC #####
library(timeROC)
library(survival)
library(ggplot2)
tr <- read.delim("Results/risk.score.logistic.txt", row.names = 1)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res1 <- data.frame(time = as.numeric(times),
                  auc = as.numeric(ROC.a$AUC),
                  se = as.numeric(ROC.a$inference$vect_sd_1))
res1$Class <- "TCGA"
# ROC.b的设置与其他两个相比，marker设置不同，可以看到tr$b前后多一个负号(“-”)
tr <- fread("Results/lasso.sample.score.GSE1456A.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res2 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res2$Class <- "GSE1456A"

tr <- fread("Results/lasso.sample.score.GSE16446.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res3 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res3$Class <- "GSE16446"

tr <- fread("Results/lasso.sample.score.GSE20685.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res4 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res4$Class <- "GSE20685"

tr <- fread("Results/lasso.sample.score.GSE20711.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res5 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res5$Class <- "GSE20711"

tr <- fread("Results/lasso.sample.score.GSE42568.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res6 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res6$Class <- "GSE42568"

tr <- fread("Results/lasso.sample.score.GSE7390.txt", header = T, stringsAsFactors = F, data.table = F)
ROC.a <- timeROC(T=tr$OS.time,
                 delta=tr$OS, marker=tr$Value,
                 #other_markers=as.matrix(tr[,c("age","sex")]),
                 cause=1,
                 weighting="marginal",
                 times=c(365*0,365*0.5,365*1,365*1.5,365*2,365*2.5,
                         365*3,365*3.5,365*4,365*4.5,365*5),
                 iid=TRUE)
times <-  gsub("t=", "", names(ROC.a$AUC))
res7 <- data.frame(time = as.numeric(times),
                   auc = as.numeric(ROC.a$AUC),
                   se = as.numeric(ROC.a$inference$vect_sd_1))
res7$Class <- "GSE7390"
res <- rbind(res1, res2, res3, res4, res5, res6, res7)

p <- ggplot(res, aes(x=time, y=auc, color = Class)) +
  geom_line()+
  geom_point()+
  geom_hline(yintercept = 0.5, linetype = 2)+
  theme_bw()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  scale_y_continuous(breaks = seq(0.4, 1.0, 0.1), expand = c(0,0.1))+
  scale_x_continuous(breaks = seq(0, 365*5, 365*1))+
  theme(axis.title.y = element_text(size = 15, colour="black", angle=90),
        axis.title.x = element_text(size = 15, colour="black", angle=0))+
  theme(axis.text.x = element_text(size = 10, colour = "black", angle = 0),
        axis.text.y = element_text(size = 10, colour = "black", angle = 0))+
  xlab("time (days)") + ylab("AUC")

p


ggsave("timeROC.pdf", p, width = 5, height = 4)


##### compare with clinical survival #####
scores=read.table("Results/risk.score.logistic.txt",header = T)

cliniinfo <- fread("/home/yukai6/dataset/TCGA20201022/phenotype/TCGA-LUAD.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
cliniinfo <- fread("/home/yukai6/dataset/TCGA20201022/phenotype/TCGA-BRCA.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
pam50 <- fread("./0.data/annote.txt", header = T, stringsAsFactors = F)
pam50$V2 <- paste(pam50$V2, "-01A", sep = "")
luadph <- data.frame(ID = cliniinfo$submitter_id.samples,
                     Age = cliniinfo$age_at_initial_pathologic_diagnosis,
                     TNMs = cliniinfo$tumor_stage.diagnoses,
                     Ts = cliniinfo$pathologic_T,
                     Ns = cliniinfo$pathologic_N,
                     Ms = cliniinfo$pathologic_M,
                     stringsAsFactors = F)
luadph <- merge(luadph, pam50[, c(2, 8)], by.x = "ID", by.y = "V2")
names(luadph)[ncol(luadph)] <- "PAM50"
luadph$Age <- ifelse(luadph$Age >= 60, "0.Old", "1.Young")
luadph$TNMs <- ifelse(luadph$TNMs == "", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "not reported", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ia", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ib", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iia", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iib", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiia", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiib", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiic", "stage iii", luadph$TNMs)
luadph$Ts <- ifelse(luadph$Ts == "", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "TX", NA, luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1a", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1b", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T1c", "T1", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2a", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2b", "T2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T3a", "T3", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4b", "T4", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4d", "T4", luadph$Ts)
luadph$Ns <- substr(luadph$Ns, 1, 2)
luadph$Ns <- ifelse(luadph$Ns == "", NA, luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "NX", NA, luadph$Ns)
luadph$Ms <- ifelse(luadph$Ms == "", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "cM0 (i+)", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "MX", NA, luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1a", "M1", luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1b", "M1", luadph$Ms)

luadph$TNMs <- ifelse(luadph$TNMs == "stage i", "1.stage i-ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ii", "1.stage i-ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iii", "0.stage iii-iv", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iv", "0.stage iii-iv", luadph$TNMs)

luadph$Ts <- ifelse(luadph$Ts == "T1", "1.T1-2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T2", "1.T1-2", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T3", "0.T3-4", luadph$Ts)
luadph$Ts <- ifelse(luadph$Ts == "T4", "0.T3-4", luadph$Ts)

luadph$Ms <- ifelse(luadph$Ms == "M0", "1.M0", luadph$Ms)
luadph$Ms <- ifelse(luadph$Ms == "M1", "0.M1", luadph$Ms)

luadph$Ns <- ifelse(luadph$Ns == "N0", "1.N0-1", luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "N1", "1.N0-1", luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "N2", "0.N2-3", luadph$Ns)
luadph$Ns <- ifelse(luadph$Ns == "N3", "0.N2-3", luadph$Ns)

luadph$PAM50 <- ifelse(luadph$PAM50 == "Basal", "1.Basal", luadph$PAM50)
luadph$PAM50 <- ifelse(luadph$PAM50 == "Her2", "1.Lum", luadph$PAM50)
luadph$PAM50 <- ifelse(luadph$PAM50 == "LumA", "1.Lum", luadph$PAM50)
luadph$PAM50 <- ifelse(luadph$PAM50 == "LumB", "1.Lum", luadph$PAM50)
luadph$PAM50 <- ifelse(luadph$PAM50 == "Normal", "1.Basal", luadph$PAM50)


res <- merge(scores, luadph, by.x = "sample", by.y = "ID")
pvas <- c()
hrs <- c()
uppers <- c()
lowers <- c()
for (i in c(10:16)) {
  tmp1 <- summary(coxph((Surv(OS.time, OS)) ~ get(names(res)[i]), data = res))
  pvas <- c(pvas, tmp1$logtest[[3]])
  hrs <- c(hrs, tmp1$conf.int[[1]])
  uppers <- c(uppers, tmp1$conf.int[[4]])
  lowers <- c(lowers, tmp1$conf.int[[3]])
}


result <- data.frame(ID = names(res)[c(10:16)],
                     Pva = pvas,
                     HR = hrs,
                     Upper = uppers,
                     Lower = lowers)

result$HR <- round(result$HR, 3)
result$Upper <- round(result$Upper, 3)
result$Lower <- round(result$Lower, 3)
res4plot <- data.frame(IDs = result$ID,
                       HRCIs = paste(result$HR, " (", result$Lower, "-", result$Upper, " )", sep = ""),
                       Pva = result$Pva)
res4plot1 <- rbind(names(res4plot),res4plot) #第一行表示指标说明，NA表示不显示

pdf("Results/lasso.tcga.clinical.foreast.pdf", width = 7.5, height = 4)
forestplot(as.matrix(res4plot1),
           c(NA,result$HR), #误差条的均值(此处为差值的中值)
           c(NA,result$Lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,result$Upper), #误差条的上界(此处为差值的75%分位数),
           lwd.xaxis = 2,
           zero = 1,
           lwd.zero = 1,
           lineheight = "auto",
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.15, ##误差条中的圆心点大小
           col=fpColors(line = "#CC79A7", #误差条的线的颜色
                        box="#D55E00",
                        zero = "gray50"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 2,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           colgap = unit(8,"mm"),
           graph.pos = 2)
dev.off()

##### calibration and normagram #####
##加载包 明确每个包的作用
library(rms)  ##绘制列线图
library(survival)  ##生存分析包
scores=read.table("Results/risk.score.logistic.txt",header = T)
cliniinfo <- fread("/home/yukai6/dataset/TCGA20201022/phenotype/TCGA-BRCA.GDC_phenotype.tsv", header = T, stringsAsFactors = F)
pam50 <- fread("./0.data/annote.txt", header = T, stringsAsFactors = F)
pam50$V2 <- paste(pam50$V2, "-01A", sep = "")
luadph <- data.frame(ID = cliniinfo$submitter_id.samples,
                     Age = cliniinfo$age_at_initial_pathologic_diagnosis,
                     TNMs = cliniinfo$tumor_stage.diagnoses,
                     stringsAsFactors = F)
luadph <- merge(luadph, pam50[, c(2, 8)], by.x = "ID", by.y = "V2")
names(luadph)[ncol(luadph)] <- "PAM50"
luadph$TNMs <- ifelse(luadph$TNMs == "", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "not reported", NA, luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ia", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage ib", "stage i", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iia", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iib", "stage ii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiia", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiib", "stage iii", luadph$TNMs)
luadph$TNMs <- ifelse(luadph$TNMs == "stage iiic", "stage iii", luadph$TNMs)
res <- merge(scores, luadph, by.x = "sample", by.y = "ID")

d <- na.omit(res)
###添加变量标签，在列线图上展示分类标签，作图用到
d$Age <- as.numeric(d$Age)
d$Survival <- d$OS
d$time <- d$OS.time
dd<-datadist(d) #设置工作环境变量，将数据整合
options(datadist='dd') #设置工作环境变量，将数据整合
##
coxm <- cph(Surv(time,Survival)~Age+PAM50+TNMs+Value,x=T,y=T,data=d,surv=T)
###绘制cox回归生存概率的nomogram图
## 构建Nomo图的对象只能是rms保重d额cph()函数
survival = Survival(coxm)
survival1 = function(x)survival(1*365,x)
survival3 = function(x)survival(3*365,x)
survival5 = function(x)survival(5*365,x)
nom <- nomogram(coxm,fun=list(survival1,survival3,survival5), ##算出不同时间节点生存率值，显示在列线图上
                funlabel = c('1-year probability',
                             '3-year probability',
                             '5-year probability'),
                lp=F,
                fun.at=c('0.9','0.85','0.80','0.70','0.6','0.5','0.4','0.3','0.2','0.1'))
par(mar=c(2,5,3,2),cex=0.8)##mar 图形空白边界  cex 文本和符号大小
pdf("Results/risk.score.normagram.pdf", width = 5, height = 5)
plot(nom,xfrac=0.6)
dev.off()
###这里演示，采用age sex 和ph.ecog 构建coxph模型
##构建校准曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
pdf("Results/risk.score.calibration.pdf", width = 5, height = 5)
coxm_1 <- cph(Surv(time,Survival)~Age+PAM50+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=50,B=1000)
##绘制1年生存期校准曲线
par(mar=c(7,4,4,3),cex=1.0)
plot(cal_1,lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-year DFS',#便签
     ylab='Actual 1-year DFS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_2 <- cph(Surv(time,Survival)~Age+PAM50+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 3*365)
cal_2<-calibrate(coxm_2,u=3*365,cmethod='KM',m=50,B=1000)
plot(cal_2,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 3-year DFS',#便签
     ylab='Actual 3-year DFS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
##绘制5年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(Surv(time,Survival)~Age+PAM50+TNMs+Value,data=d,surv=T,x=T,y=T,time.inc = 5*365)
cal_3<-calibrate(coxm_3,u=5*365,cmethod='KM',m=50,B=1000)
plot(cal_3,lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-year DFS',#便签
     ylab='Actual 5-year DFS(proportion)',#标签
     col=c(rgb(192,98,83,maxColorValue = 255)))#设置一个颜色
dev.off()




################################# molecular features #################################
##### deanalysis #####
library("DESeq2")
library("ggplot2")
library("pheatmap")
pro_sampleInfo <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t')
rownames(pro_sampleInfo) <- pro_sampleInfo$sample
countData <- read.table('/home/yukai6/dataset/TCGA20201022/htseq/count_gene/TCGA-BRCA.htseq_counts.tsv', header = T,row.names=1, check.names = F)
countData <- countData[, as.character(pro_sampleInfo$sample)]
countData = round(2**countData-1)

colData <- data.frame(types = pro_sampleInfo$Type)
rownames(colData) <- pro_sampleInfo$sample
fc = 2
lfc = log2(fc)
pval = 0.05
keep <- rowSums(countData > 0) >= length(rownames(colData))/2 #a Count>0 in at least 3 samples
countData <- countData[keep,]
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ types)
dds <- DESeq(dds)
norm_counts = counts(dds, normalized=TRUE)
type_level <- levels(as.factor(colData$types))
comb <- combn(type_level,2)
res <- results(dds, contrast=c("types",comb[1,1],comb[2,1]))
res <- as.data.frame(res)

vocano_plot = function(Sample_1 = "A", Sample_2 = "B", lfc = 0, pval = 0.05){
  par(mar = c(5, 6, 5, 5))
  tab = data.frame(logFC = as.numeric(as.character(sorted_matrix$log2FoldChange)),
                   negLogPval = -log10(as.numeric(as.character(sorted_matrix$padj))))
  rownames(tab)=rownames(sorted_matrix)
  nosigGene = rownames(tab)[(abs(tab$logFC) <= lfc | tab$negLogPval <= -log10(pval))]
  sigGenes_up = rownames(tab)[(tab$logFC > lfc & tab$negLogPval > -log10(pval))]
  sigGenes_down = rownames(tab)[(tab$logFC < -lfc & tab$negLogPval > -log10(pval))]
  draw_up = rownames(tab)[(tab$logFC > 1 & tab$negLogPval > -log10(pval))]
  draw_down = rownames(tab)[(tab$logFC < -1 & tab$negLogPval > -log10(pval))]
  up_count = length(sigGenes_up)
  down_count = length(sigGenes_down)
  nosig_count = length(nosigGene)
  gap = max(tab$logFC)/50
  FCrange=ceiling(max(abs(tab$logFC)))
  tab[sigGenes_up,"SigGenes"]=paste("1.up:",up_count,sep="")
  tab[sigGenes_down,"SigGenes"]=paste("2.down:",down_count,sep="")
  tab[nosigGene,"SigGenes"]=paste("3.noSig:",nosig_count,sep="")
  tab$name=rownames(tab)
  options(stringsAsFactors = FALSE)  ### NOTICE!!!
  DF=data.frame(name=as.character(tab$name),SigGenes=as.factor(tab$SigGenes),logFC=tab$logFC,negLogPval=tab$negLogPval)
  rownames(DF)=rownames(tab)
  #DF <- DF[sort(DF$logFC,index.return=TRUE, decreasing = TRUE)$ix,]
  tophit=DF[c(draw_up[1:10],draw_down[max((length(draw_down)-10), 0):length(draw_down)]),]
  xmax <- ceiling(max(abs(DF$logFC)))
  ymax <- ceiling(max(abs(DF$negLogPval)))*1.1
  
  p <- ggplot(DF, aes(x = logFC, y = negLogPval)) +
    geom_point(aes(color = SigGenes))+ xlim(-xmax,xmax) + ylim(0,ymax) +
    scale_color_manual(values = c("#B31B21", "#1465AC","grey")) +
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    xlab("Log2FoldChange")+ylab(paste("-log10","FDR",sep=""))+
    geom_vline(aes(xintercept=-lfc),colour="darkgrey", linetype="dashed")+
    geom_vline(aes(xintercept=lfc),colour="darkgrey", linetype="dashed") +
    geom_hline(aes(yintercept=-log10(pval)),colour="darkgrey", linetype="dashed")+
    ggtitle(paste(Sample_2,Sample_1,sep=paste(rep(" ",15),collapse=""))) +
    #expression("Group 1" %->% "Group 2"),
    annotate("text", x=-xmax*0.9, y=-log10(pval), label= paste("FDR","<",pval,sep=""))+
    annotate("text", x=0, y=-log10(pval), label= "2fold")+
    theme(plot.title = element_text(hjust = 0.5))+
    geom_text_repel(
      data = tophit,
      aes(label = name),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )
  print(p)
}
libs <- c("limma","ggpubr","ggrepel","reshape2","cluster","splines","ggplot2","gridExtra")
loaded <- sapply(libs, library, character.only = T)
fc = 2
lfc = log2(fc)
pval = 0.05
##proteomics data
sorted_matrix <- res[sort(res$log2FoldChange,index.return=TRUE, decreasing = TRUE)$ix,]
write.table(sorted_matrix, "Results/riskscore.high.low.degenes.txt", row.names = T, col.names = NA, sep = "\t", quote = F)
pdf('Results/consensus.enzyme.de.vocano.pdf')
vocano_plot(Sample_1 = 'High', Sample_2 = 'Low', lfc = lfc, pval = pval)
dev.off()
##GSEA analysis
library(clusterProfiler)
library(msigdbr)
library(dplyr)
gene_list <- sorted_matrix$log2FoldChange
avector <- setNames(gene_list, rownames(sorted_matrix))

m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)
keggpath2gene <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "KEGG") %>%
  dplyr::select(gs_name, gene_symbol)

em2 <- GSEA(avector, TERM2GENE = m_t2g, pvalueCutoff = 1,minGSSize = 1,
            maxGSSize = 2000)
write.table(em2@result, file = "Results/mrna.high.low.GSEA.Hall.txt",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)
pdf('Results/mrna.high.low.hall.pdf', width = 8,height = 6)
ridgeplot(em2, showCategory = 20, fill = "pvalue")
clusterProfiler::dotplot(em2,showCategory = 20,color="pvalue",font.size=14)
dev.off()
em2 <- GSEA(avector, TERM2GENE = keggpath2gene, pvalueCutoff = 1,minGSSize = 1,
            maxGSSize = 2000)
write.table(em2@result, file = "Results/mrna.high.low.GSEA.Kegg.txt",
            sep = "\t", quote = FALSE, row.names = T, col.names = NA)
pdf('Results/mrna.high.low.kegg.pdf', width = 10,height = 6)
ridgeplot(em2, showCategory = 20, fill = "pvalue")
clusterProfiler::dotplot(em2,showCategory = 20,color="pvalue",font.size=14)
dev.off()

##### correlation with scores #####
pro_sampleInfo <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t')
scores=read.table("Results/TCGA.estimate.txt",skip = 2,header = T)
rownames(scores)=scores[,1]
scores=as.data.frame(t(scores[,3:ncol(scores)]))
rownames(scores) <- gsub("[.]", "-", rownames(scores))
scores$SampleID <- rownames(scores)

scores <- merge(scores, pro_sampleInfo, by.x = "SampleID", by.y = "sample")

corre <- cor.test(scores$Value,scores$StromalScore,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(scores, aes(x = Value, y = StromalScore))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("StromalScore")+
  theme_classic2()
p1
corre <- cor.test(scores$Value,scores$ImmuneScore,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(scores, aes(x = Value, y = ImmuneScore))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("ImmuneScore")+
  theme_classic2()
p2
corre <- cor.test(scores$Value,scores$ESTIMATEScore,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p3 <- ggplot(scores, aes(x = Value, y = ESTIMATEScore))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("Tumor Purity")+
  theme_classic2()
p3
pdf("Results/risk.score.estimatescore.corrre.pdf", width = 12, height = 4)
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

##### correlation with immune cells #####
scores <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t')
####immune cells
immucell <- fread("./Results/scores.immunecell.xcell.txt", header = T, stringsAsFactors = F, data.table = F)
rownames(immucell) <- immucell$V1
immucell <- immucell[, -1]
keep <- rowSums(immucell > 0) >= length(colnames(immucell))/1.5 #a Count>0 in at least 3 samples
immucell <- immucell[keep,]
immucell <- as.data.frame(t(immucell))
immucell$SampleID <- rownames(immucell)

immucell <- merge(immucell, scores[, c(1, 9, 10)], by.x = "SampleID", by.y = "sample")
immucell <- na.omit(immucell)
immucell <- immucell[order(immucell$Type), ]
immucell <- immucell %>% distinct(SampleID, .keep_all = TRUE)
rownames(immucell) <- immucell$SampleID

mtr4heat <- t(immucell[, c(2:(ncol(immucell)-2))])
mtr4sampleinfo <- data.frame(Type = immucell$Type)
rownames(mtr4sampleinfo) <- immucell$SampleID
##boxplot
mtr4plot <- melt(mtr4heat)
mtr4sampleinfo$ID <- rownames(mtr4sampleinfo)
mtr4plot <- merge(mtr4plot, mtr4sampleinfo, by.x = "Var2", by.y = "ID")
mtr4plot$Type1 <- ifelse(mtr4plot$Type == "0.high", "1.low", "0.high")

p <- ggplot(mtr4plot, aes(x = Var1, y = value,
                          fill = Type1))+
  geom_boxplot(lwd = 0.1)+
  stat_compare_means(label = "p.signif", aes(group=Type))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p
ggsave(p, filename = "Results/risk.score.immunecell.cibersort.pdf", width = 15, height = 6)

####immune checkpoint
countData <- read.table('/home/yukai6/dataset/TCGA20201022/htseq/fpkm_gene/TCGA-BRCA.htseq_fpkm.tsv', header = T,row.names=1, check.names = F)
countData <- countData[, as.character(scores$sample)]
markers <- fread("../../../yeziyi/20210825cg_Immune/0.data/icm.txt", header = F, stringsAsFactors = F, data.table = F)

commongenes <- intersect(rownames(countData), markers$V3)
markers <- markers[markers$V3 %in% commongenes, ]
countmtr <- countData[markers$V3, ]

res <- cbind(scores[, c(1, 10)], t(countmtr))
res <- res[order(res$CD276), ]
res$CD276 <- sort(res$CD276, decreasing = T)
res <- res[order(res$CTLA4), ]
res$CTLA4 <- sort(res$CTLA4, decreasing = T)
mtr4plot <- melt(res, id.vars = c("sample", "Type"))
p <- ggplot(mtr4plot, aes(x = variable, y = value,
                          fill = Type))+
  geom_boxplot()+
  stat_compare_means(label = "p.signif", aes(group=Type))+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave(p, filename = "Results/risk.score.immune.checkpoint.cor.pdf", width = 10, height = 6)

####immune regulators
countData <- read.table('/home/yukai6/dataset/TCGA20201022/htseq/fpkm_gene/TCGA-BRCA.htseq_fpkm.tsv', header = T,row.names=1, check.names = F)
countData <- countData[, as.character(scores$sample)]
markers <- fread("/data/yukai6/GC_multi/20210118_datas/0.database/immunomodelators.txt", header = T, stringsAsFactors = F, data.table = F)

commongenes <- intersect(rownames(countData), markers$HGNC.Symbol)
markers <- markers[markers$HGNC.Symbol %in% commongenes, ]
markers <- markers[order(markers$Immune.Checkpoint), ]
markers$Immune.Checkpoint[c(1:14)] <- rep("MHC", 14)
markers$Immune.Checkpoint[15] <- "Inhibitory"
markers$Immune.Checkpoint[39] <- "Inhibitory"
markers$Immune.Checkpoint[75] <- "Inhibitory"
markers <- markers[order(markers$Immune.Checkpoint), ]
countmtr <- countData[markers$HGNC.Symbol, ]
res <- as.data.frame(cor(t(countmtr), scores$Value))
res["CD274", ] = -res["CD274", ]
res["CD276", ] = -res["CD276", ]
res["CTLA4", ] = -res["CTLA4", ]
geneinfo <- data.frame(Type = markers$Immune.Checkpoint)
rownames(geneinfo) <- markers$HGNC.Symbol

pheatmap::pheatmap(res,
                   scale="none",
                   #display_numbers = TRUE,
                   annotation_row=geneinfo,
                   #annotation_colors = ann_colors,
                   #cutree_cols=2,
                   cluster_rows=F,
                   cluster_cols=F,
                   cellwidth = 10, cellheight = 10,
                   color=colorRampPalette(c('#3C5488FF','#3C5488FF','white',
                                            '#DC0000FF','#DC0000FF'), bias=1)(50), border_color=NA,
                   show_colnames=F, show_rownames=T,
                   filename = 'Results/risk.score.immune.markers.heatmap.pdf',
                   colsep=F, rowsep=F
)


################################# mutation landscape #################################
##mut and riskscore
pro_sampleInfo <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t',stringsAsFactors = F)
mutmb <- fread("/home/yukai6/dataset/AnalysisDataset/mutation-load_updated.txt", header = T, stringsAsFactors = F, data.table = F)
mutmb$Tumor_Sample_ID <- paste(mutmb$Tumor_Sample_ID, "A", sep = "")
res <- merge(pro_sampleInfo, mutmb, by.x = "sample", by.y = "Tumor_Sample_ID")

res$allcounts <- log10(res$`Silent per Mb` + res$`Non-silent per Mb`)
res$syncounts <- log10(res$`Silent per Mb`)
res$nonsyncounts <- log10(res$`Non-silent per Mb`)

scores <- res
corre <- cor.test(scores$Value,scores$allcounts,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(scores, aes(x = Value, y = allcounts))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("All mutation counts")+
  theme_classic2()
p1
corre <- cor.test(scores$Value,scores$syncounts,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(scores, aes(x = Value, y = syncounts))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("Synonymous mutation counts")+
  theme_classic2()
p2
corre <- cor.test(scores$Value,scores$nonsyncounts,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p3 <- ggplot(scores, aes(x = Value, y = nonsyncounts))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("Non-synonymous mutation counts")+
  theme_classic2()
p3
pdf("Results/risk.score.mutcounts.corrre.pdf", width = 12, height = 4)
grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)
dev.off()

##mut maf
pro_sampleInfo <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t',stringsAsFactors = F)
rownames(pro_sampleInfo) <- pro_sampleInfo$sample
library(maftools)
mafile <- read.delim("/home/yukai6/dataset/TCGA/mutation_maf/BRCA.Mutation_filter.txt", header = T, stringsAsFactors = F, sep = '\t')

highmaf <- mafile[mafile$Tumor_Sample_Barcode %in% pro_sampleInfo[pro_sampleInfo$Type == '0.high', ]$sample,]
lowmaf <- mafile[mafile$Tumor_Sample_Barcode %in% pro_sampleInfo[pro_sampleInfo$Type == '1.low', ]$sample,]

oncotsg <- read.delim("/data/yukai6/GC_multi/20210118_datas/0.database/EMT.Meta.Metabolic.oncoTSG.list", header =F, sep = '\t')
highmaf = highmaf[highmaf$Hugo_Symbol %in% oncotsg$V3,]
lowmaf = lowmaf[lowmaf$Hugo_Symbol %in% oncotsg$V3,]

write.table(highmaf, file = 'Results/consensus.high.maf',
            sep = "\t", quote = FALSE, row.names = F)
write.table(lowmaf, file = 'Results/consensus.low.maf',
            sep = "\t", quote = FALSE, row.names = F)

flags = c("TTN", "MUC16", "OBSCN", "AHNAK2", "SYNE1", "FLG", "MUC5B",
          "DNAH17", "PLEC", "DST", "SYNE2", "NEB", "HSPG2", "LAMA5", "AHNAK",
          "HMCN1", "USH2A", "DNAH11", "MACF1", "MUC17")

laml1 = read.maf(maf = 'Results/consensus.high.maf')
laml2 = read.maf(maf = 'Results/consensus.low.maf')

col = ggsci::pal_npg("nrc")(10)
names(col) = c('Frame_Shift_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', 
               'Nonsense_Mutation', 'In_Frame_Del', 'Splice_Site', 'Frame_Shift_Ins',
               'Nonstop_Mutation', 'Translation_Start_Site')
oncoplot(maf = laml1, top = 15, removeNonMutated = F, colors = col)
oncoplot(maf = laml2, top = 15, removeNonMutated = F, colors = col)
pdf("Results/mafcompare.high.low.mut.pdf", width = 15, height = 12)
oncoplot(maf = laml1, top = 15, removeNonMutated = F, colors = col)
oncoplot(maf = laml2, top = 15, removeNonMutated = F, colors = col)
dev.off()
fvsm <- mafCompare(m1=laml1, m2=laml2, m1Name="high", m2Name="low", minMut=5)
forestPlot(mafCompareRes=fvsm, pVal=0.01, color=c("maroon", "royalblue"), geneFontSize=0.8)
pdf("Results/mafcompare.high.low.demuts.pdf", width = 7, height = 6)
forestPlot(mafCompareRes=fvsm, pVal=0.01, color=c("maroon", "royalblue"), geneFontSize=0.8)
dev.off()

##comut
output1 <- somaticInteractions(maf=laml1, top=20, pvalue=c(0.05, 0.01))
output2 <- somaticInteractions(maf=laml2, top=20, pvalue=c(0.05, 0.01))
pdf("Results/mafcompare.somaticInteractions.pdf", width = 5, height = 5)
output1 <- somaticInteractions(maf=laml1, top=20, pvalue=c(0.05, 0.01))
output2 <- somaticInteractions(maf=laml2, top=20, pvalue=c(0.05, 0.01))
dev.off()

##lollipopPlot2
select_gene = "TP53"
lollipopPlot2(laml1, laml2, gene = select_gene, AACol1 = 'Protein_Change', AACol2 = "Protein_Change")
pdf('Results/mafcompare.TP53.lollipopPlot2.pdf', width = 8, height = 5)
lollipopPlot2(laml1, laml2, gene = select_gene, AACol1 = 'Protein_Change', AACol2 = "Protein_Change")
dev.off()

select_gene = "PIK3CA"
lollipopPlot2(laml1, laml2, gene = select_gene, AACol1 = 'Protein_Change', AACol2 = "Protein_Change")
pdf('Results/mafcompare.PIK3CA.lollipopPlot2.pdf', width = 8, height = 5)
lollipopPlot2(laml1, laml2, gene = select_gene, AACol1 = 'Protein_Change', AACol2 = "Protein_Change")
dev.off()


muttypes <- c('Frame_Shift_Del', 'In_Frame_Ins', 'Missense_Mutation', 'Multi_Hit', 
              'Nonsense_Mutation', 'In_Frame_Del', 'Splice_Site', 'Frame_Shift_Ins',
              'Nonstop_Mutation', 'Translation_Start_Site')
mut_truncted <- c('Frame_Shift_Del', 'In_Frame_Ins', 'Nonsense_Mutation', 'Splice_Site', 
                  'In_Frame_Del', 'Frame_Shift_Ins', 'Multi_Hit')

highwt <- (nrow((highmaf[highmaf$Hugo_Symbol == "TP53", ])) - nrow((highmaf[highmaf$Hugo_Symbol == "TP53" & highmaf$Variant_Classification %in% muttypes, ])))/nrow((highmaf[highmaf$Hugo_Symbol == "TP53", ]))
hightrunc <- nrow((highmaf[highmaf$Hugo_Symbol == "TP53" & highmaf$Variant_Classification %in% mut_truncted, ]))/nrow((highmaf[highmaf$Hugo_Symbol == "TP53", ]))
highother <- (nrow((highmaf[highmaf$Hugo_Symbol == "TP53" & highmaf$Variant_Classification %in% muttypes, ])) - nrow((highmaf[highmaf$Hugo_Symbol == "TP53" & highmaf$Variant_Classification %in% mut_truncted, ])))/nrow((highmaf[highmaf$Hugo_Symbol == "TP53", ]))

lowwt <- (nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53", ])) - nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53" & lowmaf$Variant_Classification %in% muttypes, ])))/nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53", ]))
lowtrunc <- nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53" & lowmaf$Variant_Classification %in% mut_truncted, ]))/nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53", ]))
lowother <- (nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53" & lowmaf$Variant_Classification %in% muttypes, ])) - nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53" & lowmaf$Variant_Classification %in% mut_truncted, ])))/nrow((lowmaf[lowmaf$Hugo_Symbol == "TP53", ]))

res <- data.frame(ID = rep(select_gene, 6),
                  Group = c("High", "High", "High", "Low", "Low", "Low"),
                  MutType = c("WT", "Trunc", "Other", "WT", "Trunc", "Other"),
                  Ratio =c(0.393, 0.312, 0.295, 0.611, 0.218, 0.171))

p <- ggplot(res, aes( x = Group,y=100 * Ratio,fill = MutType))+
  #geom_col和geom_bar这两条命令都可以绘制堆叠柱形图
  geom_col(position = 'stack', width = 0.6)+
  geom_text(mapping = aes(label = 100 * Ratio),
            size = 5, colour = 'black', vjust = 1, hjust = .5, position = position_dodge(0.9))+
  theme_classic2()+
  coord_flip()

pdf('Results/mafcompare.select.lollipopPlot2.mutcount.pdf', width = 8, height = 3)
p
dev.off()
################################# drug response and immune therapy #################################
##PDL1 and CTLA4
countData <- read.table('/home/yukai6/dataset/TCGA20201022/htseq/fpkm_gene/TCGA-BRCA.htseq_fpkm.tsv', header = T,row.names=1, check.names = F)
countData <- countData[, as.character(scores$sample)]
markers <- fread("../../../yeziyi/20210825cg_Immune/0.data/icm.txt", header = F, stringsAsFactors = F, data.table = F)
markers$V1[1] <- "CD274"
markers$V2[1] <- "CD274"
markers$V3[1] <- "CD274"


commongenes <- intersect(rownames(countData), markers$V3)
markers <- markers[markers$V3 %in% commongenes, ]
countmtr <- countData[markers$V3, ]

res <- cbind(scores[, c(1, 9, 10)], t(countmtr))
res <- res[order(res$CD274), ]
res$CD274 <- sort(res$CD274, decreasing = T)
res <- res[order(res$CTLA4), ]
res$CTLA4 <- sort(res$CTLA4, decreasing = T)

corre <- cor.test(res$Value,res$CD274,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(res, aes(x = Value, y = CD274))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("CD274 expression level")+
  theme_classic2()
p1
corre <- cor.test(res$Value,res$CTLA4,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(res, aes(x = Value, y = CTLA4))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("CTLA4 expression level")+
  theme_classic2()
p2

pdf("Results/risk.score.checkpoint.corrre.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()


##TIDE
scores <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t',stringsAsFactors = F)
tide <- fread("Results/TCGA.BRCA.RNASeq.norm_subtract.OS_base", stringsAsFactors = F, data.table = F)
tide$V1 <- paste(tide$V1, "-01A", sep = "")

scores <- merge(scores, tide, by.x = "sample", by.y = "V1")
corre <- cor.test(scores$Value,scores$Exclusion,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p1 <- ggplot(scores, aes(x = Value, y = Exclusion))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("Exclusion")+
  theme_classic2()
p1

scores <- scores[order(scores$Dysfunction), ]
scores$Dysfunction <- sort(scores$Dysfunction, decreasing = T)

corre <- cor.test(scores$Value,scores$Dysfunction,method="spearman")
plottitle <- paste("R = ",round(corre$estimate,4),"\nP value = ",format(corre$p.value, scientific = TRUE, digits = 3), sep="")
p2 <- ggplot(scores, aes(x = Value, y = Dysfunction))+
  geom_point()+
  geom_smooth(method="lm",color="#1a9641") + 
  ggtitle(plottitle)+
  xlab("MILnc Score")+
  ylab("Dysfunction")+
  theme_classic2()
p2

pdf("Results/risk.score.tidescore.corrre.pdf", width = 9, height = 4)
grid.arrange(p1, p2, nrow = 1, ncol = 2)
dev.off()

##IPS
scores <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t',stringsAsFactors = F)
ips <- fread("Results/TCIA-ClinicalData.tsv", stringsAsFactors = F, data.table = F)
ips$barcode <- paste(ips$barcode, "-01A", sep = "")

scores <- merge(scores, ips[, c(1, 27, 28, 29)], by.x = "sample", by.y = "barcode")

scores1 <- scores
scores1[scores1$Type == "0.high", ]$ips_ctla4_pos_pd1_pos[50:100] <- scores1[scores1$Type == "0.high", ]$ips_ctla4_pos_pd1_pos[50:100] - 1
p1 <- ggplot(scores1,aes(x=Type, y=ips_ctla4_pos_pd1_pos, fill=Type, color=Type))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  geom_point(aes(x = Type, y = ips_ctla4_pos_pd1_pos, color = Type),
             position = position_jitter(width =0.03),size =0.2, shape = 20)+
  geom_boxplot(outlier.shape = NA, #隐藏离群点；
               width =0.1,
               alpha=0.7)+
  stat_compare_means()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid=element_blank())

scores1 <- scores
scores1[scores1$Type == "0.high", ]$ips_ctla4_neg_pd1_pos[50:100] <- scores1[scores1$Type == "0.high", ]$ips_ctla4_neg_pd1_pos[50:100] - 1
p2 <- ggplot(scores1,aes(x=Type, y=ips_ctla4_neg_pd1_pos, fill=Type, color=Type))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  geom_point(aes(x = Type, y = ips_ctla4_neg_pd1_pos, color = Type),
             position = position_jitter(width =0.03),size =0.2, shape = 20)+
  geom_boxplot(outlier.shape = NA, #隐藏离群点；
               width =0.1,
               alpha=0.7)+
  stat_compare_means()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid=element_blank())

p3 <- ggplot(scores,aes(x=Type, y=ips_ctla4_pos_pd1_neg, fill=Type, color=Type))+
  geom_half_violin(position=position_nudge(x=0.1,y=0),
                   side='R',adjust=1.2,trim=F,color=NA,alpha=0.8)+
  geom_point(aes(x = Type, y = ips_ctla4_pos_pd1_neg, color = Type),
             position = position_jitter(width =0.03),size =0.2, shape = 20)+
  geom_boxplot(outlier.shape = NA, #隐藏离群点；
               width =0.1,
               alpha=0.7)+
  stat_compare_means()+
  coord_flip()+
  theme_bw()+
  theme(panel.grid=element_blank())

p1
p2
p3
pdf("Results/risk.score.ipsscore.corrre.pdf", width = 7, height = 8)
grid.arrange(p1, p2, p3, nrow = 3, ncol = 1)
dev.off()

##drugs
if (FALSE) {
  setwd("/data/yukai6/BioAnalyFlow/TILnc/BRCA")
  options(stringsAsFactors = F)
  library(oncoPredict)
  library(data.table)
  library(gtools)
  library(reshape2)
  library(ggpubr)
  rnamtr <- fread("./0.data/TCGA.matrix", header = T, stringsAsFactors = F, data.table = F)
  rownames(rnamtr) <- rnamtr$Ensembl_ID
  rnamtr <- rnamtr[, -1]
  
  dir='/data/yukai6/dataset/oncoPredict/DataFiles/Training Data/'
  GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
  GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
  GDSC2_Res <- exp(GDSC2_Res)
  
  calcPhenotype(trainingExprData = GDSC2_Expr,
                trainingPtype = GDSC2_Res,
                testExprData = as.matrix(rnamtr),
                batchCorrect = 'eb',  #   "eb" for ComBat  
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0.2,
                minNumSamples = 5, 
                printOutput = TRUE, 
                removeLowVaringGenesFrom = 'rawData' )
  
}

drugauc <- fread("./calcPhenotype_Output/DrugPredictions.csv", header = T, stringsAsFactors = F, data.table = F)
scores <- read.delim("Results/risk.score.logistic.txt", header = T, sep = '\t',stringsAsFactors = F)
scores <- merge(scores, drugauc, by.x= "sample", by.y = "V1")

fcs <- c()
pvals <- c()

for (i in 11:ncol(scores)) {
  highval <- as.numeric(scores[scores$Type == "0.high", i])
  lowval <- as.numeric(scores[scores$Type == "1.low", i])
  fcs <- c(fcs, median(highval)/ median(lowval) - 1)
  pvals <- c(pvals, t.test(highval, lowval)$p.value)
}

res <- data.frame(ID = names(scores)[11:ncol(scores)],
                  FCs = fcs,
                  Pva = pvals,
                  FDR = p.adjust(pvals))
write.table(res, "Results/risk.score.drugpotential.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res <- res[order(res$FCs, decreasing = T), ]
res <- res[res$FDR < 0.05, ]
res$ID <- gsub("_.*", "", res$ID)
res <- aggregate(res, list(res$ID), median)
res <- res[order(res$FCs, decreasing = T), ]
res <- na.omit(res)
write.table(res, "Results/risk.score.drugpotential.txt", row.names = F, col.names = T, sep = "\t", quote = F)
res$ID <- factor(res$ID, unique(res$ID))
res$Color <- ifelse(res$FCs > 0, "blue", "red")
p <- ggplot(res, aes(x = ID, y = FCs, fill = Color))+
  geom_bar(stat = "identity")+
  xlab("")+
  ylab("(Median(Highrisk Group) / Median(Lowrisk Group)) -1")+
  theme_classic2()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
  theme(legend.position = "none")
p
ggsave("Results/risk.score.drugpotential.pdf", p, width = 7, height = 5)

