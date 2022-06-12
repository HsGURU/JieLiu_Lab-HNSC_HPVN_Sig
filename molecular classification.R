length(surG)
head(surG)
hpvNexp <- data.table::fread("../sup/TCGAHNSC/TNSC_SYMBOL_HPVneg.txt", 
                             data.table = FALSE) %>% 
  column_to_rownames(var = "V1")
d <- hpvNexp
rownames(d) <- str_replace_all(rownames(d), "-", ".")
setdiff(surG, rownames(d))
d <- d[surG,]
d = sweep(d,1, apply(d,1,median,na.rm=T))
d[1:5,1:5]
dim(d)
#setting output
title="../sup/TCGAHNSC/ConsensusClusterPlus/" 
results = ConsensusClusterPlus(as.matrix(d),
                               maxK=6, 
                               reps=1000,
                               pItem=0.8,
                               pFeature=1,
                               distance="pearson",
                               clusterAlg="hc",
                               #writeTable = T,
                               title=title,
                               seed=1262118388.71279,
                               plot="png")#或pdf
icl = calcICL(results,title=title,plot="png")
# take results
ablist <- lapply(2:6, function(z){
 aa <- results[[z]]$consensusClass %>% as.data.frame()
 colnames(aa)[1] <- paste0("abs", z)
 aa <- rownames_to_column(aa, var = "ID")
 aa[, 2] <- paste0("C", aa[, 2])
 aa
 })
abres <- base::Reduce(function(x, y)merge(x, y, by = "ID"), ablist)
write.csv(abres,"../sup/TCGAHNSC/CCPres.csv", row.names = FALSE)
