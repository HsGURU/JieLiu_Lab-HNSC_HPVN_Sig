# genes from WGCNA
inputgenedata <- data.frame(Symbol = names(net$colors), module = net$colors, WGCNA_module = mergedColors)
write.csv(inputgenedata, "../sup/WGCNA/input8000module.csv", row.names = F)
color = mycols[1:2]
library(survival)
library(survminer)
exp <- data.table::fread("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg_coding.txt",
                         data.table = FALSE) %>% column_to_rownames(var = "V1")
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
moduleG <- read.table("../sup/WGCNA/input8000module.csv", sep = ",", header = T)
whichmodule <- c("black", "pink",  "green", "blue", "brown")
head(moduleG)
moduleG <- moduleG %>% dplyr::filter(WGCNA_module %in% whichmodule) %>%
  dplyr::pull(Symbol)
inputG <- intersect(surG, moduleG)
length(inputG)#506
inputdata <- exp[inputG, ]  #506 302
inputdata <- scale(t(inputdata), scale = T, center = T)
trainS <- readRDS("../sup/ML/trainS_8073.rds")
testS <- readRDS("../sup/ML/testS_8073.rds")
trainS <- sample(rownames(inputdata), size = round( 0.7 * nrow(inputdata), 0),
                 replace = FALSE)
testS <- setdiff(rownames(inputdata), trainS)
traindat <- inputdata[trainS, ]
testdat <- inputdata[testS, ]
pdat <- read.table("../sup/TCGAHNSC/pd_clean_hpvN.txt", header = T, sep = "\t")
pdat <- column_to_rownames(pdat, var = "sampleID")
trainpd <- pdat[trainS, ]
testpd <- pdat[testS, ]
vfit = cv.glmnet(traindat, Surv(trainpd$OS.time,trainpd$OS), 
                  family = "cox",
                  nfold=10)
plot(cvfit)
lambda.min <- cvfit$lambda.min 
myCoefs <- coef(model, s="lambda.min") %>% as.matrix()
lasso_fea <- myCoefs[myCoefs[,1]!= 0, ,drop = FALSE] # 提出系数非0的变量
lasso_fea
length(lasso_fea)
saveRDS(lasso_fea,"../sup/ML/lasso_feature_8073.rds")
