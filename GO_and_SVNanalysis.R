##  DEGs
exp <- read.table("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg_coding.txt", header = T, sep = "\t")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
clusterdat <- read.csv("easy_input_group.csv",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
clusterdat <- filter(clusterdat, cluster %in% c("C1","C2"))
pdat <- filter(pdat,sampleID %in% rownames(clusterdat))
exp <- read.csv("../data/tcga_hnsc_mapped_expr.csv",header = T,row.names = 1)
exp <- exp[,pdat$sampleID] 
library(DESeq2)
#colData <- data.frame(row.names = rownames(clusterdat),cluter=clusterdat[,1])
dds <- DESeqDataSetFromMatrix(exp, clusterdat, design = ~cluster)
dds <- DESeq(dds)
res <- results(dds,contrast = c("cluster","C2","C1")) #C2相比C1
res <- res[order(res$pvalue),]
res <- subset(res, padj < 0.05 & abs(log2FoldChange) > 1)
head(res)
DEG <- as.data.frame(na.omit(res))
DEG <- DEG[,c(2,5)]
colnames(DEG)=c('log2FoldChange','pvalue')  
library(ggplot2)
#plot(DEG$log2FoldChange,-log2(DEG$pvalue))
logFC_cutoff <- with(DEG,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)) )
logFC_cutoff <- round(logFC_cutoff, 2)
DEG$change = as.factor(ifelse(DEG$pvalue < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
                              ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'STABLE'))
table(DEG$change)
library(EnhancedVolcano)
EnhancedVolcano(DEG,
                lab = rownames(DEG),
                x = 'log2FoldChange',
                y = 'pvalue')
save(DEG,file = "DEG_c1c2.Rdata")

## Go analysis
exp <- data.table::fread("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg_coding.txt",
                         data.table = FALSE) %>% column_to_rownames(var = "V1")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
# 
moduleG <- read.table("../sup/WGCNA/input8000module.csv", sep = ",", header = T)
whichmodule <- c("black", "pink", "greenyellow","cyan", "green", "salmon")
moduleG <- moduleG %>% dplyr::filter(WGCNA_module %in% whichmodule) 
#
write.csv(moduleG,"moduleG.csv")
moduleG <- read.csv("moduleG.csv",header = T,row.names = 1)
moduleG <- split(moduleG$Symbol, moduleG$type)
modgo.bp <- lapply(moduleG, Myenrich, category = "gobp", geneid = "SYMBOL")
c1bar <- barplot(modgo.bp$C1,showCategory =10)+ggtitle("GO biological processes enrichemnt in C1")
ggsave("../pic/gobp_C1.pdf", width = 8, height = 6)
c2bar <- barplot(modgo.bp$C2,showCategory =10)+ggtitle("GO biological processes enrichemnt in C2")
ggsave("../pic/gobp_C2.pdf", width = 8, height = 6)
library(patchwork)
p <- c1bar/c2bar
ggsave("../pic/gobp_C1C2.pdf", width = 8, height = 8)

##  SNV analysis
library(TCGAbiolinks)
library(maftools)
library(dplyr)
mut <- GDCquery_Maf(tumor = "HNSC", pipelines = "mutect2")
#clinicaldata
clinical <- GDCquery(project = "TCGA-HNSC", 
                  data.category = "Clinical", 
                  file.type = "xml")

GDCdownload(clinical)

cliquery <- GDCprepare_clinic(clinical, clinical.info = "patient")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
maf <- mutate(maf,patient_ID = str_sub(mut$Tumor_Sample_Barcode,end=16))
maf_hpvn$Tumor_Sample_Barcode <- str_replace_all(maf_hpvn$Tumor_Sample_Barcode, "-", ".")
aa <- unique(intersect(pdat$sampleID,maf_hpvn$patient_ID))#299
aa <- intersect(pdat$sampleID,snvdat$Sample_ID)
#299 patient HPV negative
maf_hpvn <- filter(maf,patient_ID %in% aa)
pdat_maf <- filter(pdat,sampleID %in% aa)
colnames(pdat_maf)[2] <- "Tumor_Sample_Barcode"
#
pd_cluster <- pd_cluster %>%
  dplyr::mutate(clusterRank  = factor(pd_clustercluster, 
                                      levels = c("C1","C2","C3")) ) %>%
  arrange(desc(cluster)) 
#
pd_cluster <- mutate(pd_cluster,cluster=factor(pd_cluster$cluster,levels = c("C1","C2","C3")))
pdf("../pic/MAF/oncoplotLIST_cluster.pdf",width = 8,height = 6)
oncoplot(maf = maf_pic,
         #genes = genelist,
         colors = col,
         annotationColor = annocolors,#给临床信息配色
         clinicalFeatures="cluster",
         sortByAnnotation = T,
         draw_titv = F)
dev.off()
#
mut <- readRDS("../data/tcga_maf.rds")
maf <- readRDS("../data/tcga_maf_hpvn.rds")
maf <- as.data.frame(maf)
jco <- c("#2874C5","#EABF00","#C6524A","#868686")
pd_cluster <- read.csv("../pd_cluster.csv",header = T ,row.names = 1)
maf_hpvn <- readRDS("../data/tcga_maf_hpvn.rds")
maf_pic <- read.maf(maf = maf_hpvn, clinicalData = pd_cluster, isTCGA = T)
maf_pic@clinical.data
#vaf
maf_pic@data$VAF <- maf_pic@data$t_alt_count/maf_pic@data$t_depth
head(maf_pic@data$Tumor_Sample_Barcode)
aa <- unique(maf_pic@data$Tumor_Sample_Barcode)
#MATH
res <- inferHeterogeneity(maf = maf_pic, tsb = aa,vafCol = 'VAF')#299patients
print(res$clusterMeans)
#pd_cluster
pdC1 <- subset(pd_cluster, cluster=="C1")$Tumor_Sample_Barcode
pdC2 <- subset(pd_cluster, cluster=="C2")$Tumor_Sample_Barcode
#maf_pic
mafC1 <- subsetMaf(maf=maf_pic, tsb=pdC1, isTCGA=TRUE)
mafC2 <- subsetMaf(maf=maf_pic, tsb=pdC2, isTCGA=TRUE)
C1vsC2 <- mafCompare(m1=mafC1, m2=mafC2, m1Name="C1", m2Name="C2")
clin_enrich <- clinicalEnrichment(maf=maf_pic, clinicalFeature="cluster")
plotEnrichmentResults(enrich_res=clin_enrich, pVal=0.01)
write.table(clin_enrich$pairwise_comparision, file="clin_enrich_pair.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(clin_enrich$groupwise_comparision, file="clin_enrich_group.tsv", quote=FALSE, row.names=FALSE, sep="\t")
#MATH score计算
resC1 <- inferHeterogeneity(maf = mafC1, tsb = 'pdC1', vafCol = 'VAF')
resC2 <- inferHeterogeneity(maf = mafC2, tsb = 'pdC2', vafCol = 'VAF')
coOncoplot(m1=mafC1, m2=mafC2, m1Name="C1", m2Name="C2")

