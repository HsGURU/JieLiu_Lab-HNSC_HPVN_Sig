exp <- data.table::fread("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg.txt",
                         data.table = FALSE) %>% 
  column_to_rownames(var = "V1")

## coding gene
probe22 <- rtracklayer::import("E:/Database/GTEx_Xena/gencode.v22.annotation.gtf") %>% as.data.frame()
probe22 <- dplyr::select(probe22, gene_id, gene_type, gene_name) %>% 
  distinct()
exp <- merge(unique(probe22[, c("gene_name", "gene_type")]), 
             exp, by.x = 1, by.y = 0)
exp <- dplyr::filter(exp, gene_type == "protein_coding") %>%
  dplyr::select(-gene_type) %>% column_to_rownames(var = "gene_name")

##  WGCNA
exp <- read.table("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg_coding.txt", header = T, sep = "\t")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
identical(colnames(exp), pdat$sampleID)
wgcnaexpr <- t(exp[order(apply(exp,1, sd), decreasing = T)[1: 8000],])
dim(wgcnaexpr)
enableWGCNAThreads()
gsg = goodSamplesGenes(wgcnaexpr, verbose = 3)
gsg$allOK
sampleTree = hclust(dist(wgcnaexpr), method = "ward.D2")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
cluster <- read.csv("../sup/TCGAHNSC/CCPres.csv", header = T)
identical(cluster$ID, pdat$sampleID)
pdat$labels <- cluster$abs3
datTraits <- pdat[, c("sampleID", "clinical_stage", "labels")]
traitColors <- datTraits[, -1]
# running WGCNA 
datExpr <- wgcnaexpr
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
pdf("../pic/WGCNA/sft.pdf", width = 6, height = 4)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
#  Picture
oftPower = sft$powerEstimate;
net = blockwiseModules(
  datExpr,
  power = sft$powerEstimate,
  maxBlockSize = 10000,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-FPKM-TOM",
  verbose = 3
)
table(net$colors)
# Hierarchical cluster dendrogram
mergedColors = labels2colors(net$colors)
table(mergedColors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datTraits$class <- factor(datTraits$labels, levels = c("C1", "C2", "C3"))
design = model.matrix(~0+ datTraits$class)
colnames(design) = levels(datTraits$class)
  
