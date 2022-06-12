library(data.table)
library(ggplot2)
library(ggsci)
library(tidyverse)
library(ggsci)
library(tableone)
library(ConsensusClusterPlus)
library(WGCNA)
library(IOBR)
library(estimate)
library(ggpubr)
library(pheatmap)
options(stringsAsFactors = FALSE)

##select HPV negative patients
clin <- data.table::fread("../data/clinical_PANCAN_patient_with_followup.tsv", data.table = FALSE)
clin <- clin[clin$acronym == "HNSC", c("bcr_patient_barcode","hpv_status_by_ish_testing", "hpv_status_by_p16_testing")]
clin$bcr_patient_barcode <- str_replace_all(clin$bcr_patient_barcode, "-", ".")
viral <- data.table::fread("../data/viral.tsv",
                           header = T, sep = "\t")
viral <- viral %>% dplyr::filter(Study == "HNSC")
viral$SampleBarcode <- str_replace_all(viral$SampleBarcode, "-", ".")
pd <- read.table("../data/pd_adj.csv", header = T, sep = ",")
pd <- pd[, -1]
pd <- dplyr::mutate(pd, patientID = str_sub(sampleID, 1, 12))
pd <- left_join(pd, viral[, c("SampleBarcode", "HPV")], 
                by = c("sampleID" = "SampleBarcode"))
pd <- left_join(pd, clin, by = c("patientID" = "bcr_patient_barcode"))
pd$hpv_by_BBT <- ifelse(pd$HPV > 10, "Positive", "Negative")


##prognistic related genes
hpvNexp <- data.table::fread("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg.txt", 
                             data.table = FALSE) %>% 
  column_to_rownames(var = "V1")
hpvNexp <- t(hpvNexp) %>% as.data.frame()
hpvNpd <- data.table::fread("../sup/TCGAHNSC/pd_clean_hpvN.txt", data.table = FALSE)
d <- data.frame(hpvNpd[, c("OS", "OS.time", "DSS", "DSS.time")], hpvNexp)
genes <- setdiff(colnames(d), c("OS", "OS.time", "DSS", "DSS.time"))
osG <- lapply(genes, Unicox, time = "OS.time", event = "OS", dat = d)
osG <- plyr::ldply(osG, "rbind") 
osG <- dplyr::mutate(osG, Type = "OS")
dssG <- lapply(genes, Unicox, time = "DSS.time", event = "DSS", dat = d)
dssG <- plyr::ldply(dssG, "rbind") %>% dplyr::mutate(.,  Type = "RFS")
aa <- dplyr::filter(osG, P < pvalue) %>% pull(gene)
bb <- dplyr::filter(dssG, P < pvalue) %>% pull(gene)
surG <- unique(c(aa, bb))
