rm(list = ls())

setwd("<dir_name>")

library(edgeR)
# coin's wilcoxon test can handle ties well
library(coin)

rawdata <- read.delim("<count_file>.txt", sep = "\t")

y <- DGEList(counts=rawdata[,3:ncol(rawdata)], genes=rawdata[,1:2])


# we remove genes that do NOT have
# at least 1 CPM in at least 5 samples (out of 87 in total), stage 4,  underweight
#keep <- filterByExpr(y)
keep <- rowSums(cpm(y)>1) >= 5
y <- y[keep, , keep.lib.sizes=FALSE]


#get normalized counts
y <- calcNormFactors(y, method="TMM")
#y$samples
normcounts<-cpm(y, normalized.lib.sizes = T)
rownames(normcounts) <- y$genes$Geneid

normcounts = as.data.frame(normcounts)

mysamples = read.csv("<sample_sheet.csv>")

rownames(mysamples) = mysamples$fastq_name

sel_samples = mysamples[colnames(normcounts),]


#############

CXCL10 = log2(normcounts['ENSG00000169245.6',]+1)
CCL5 = log2(normcounts['ENSG00000271503.7',]+1)
CXCL9 = log2(normcounts['ENSG00000138755.6',] +1)
CXCL11 = log2(normcounts['ENSG00000169248.13',] +1)
CTNNB1 = log2(normcounts['ENSG00000168036.20',] + 1)
GREM1 = log2(normcounts['ENSG00000166923.13',] + 1)




#data = data.frame(Cxcl10_log2normCounts = t(cxcl10), C1orf43_log2normCounts = t(clorf43), Rab7a_log2normCounts = t(rab7a), Sik1_log2normCounts = t(sik1), Gatd3_log2normCounts = t(gatd3), Npipb15_log2normCounts = t(npipb15))
data = data.frame(CXCL10_log2normCounts = CXCL10, CCL5_log2normCounts = CCL5, CXCL9_log2normCounts = CXCL9, CXCL11_log2normCounts = CXCL11, CTNNB1_log2normCounts = CTNNB1, GREM1_log2normCounts = GREM1)

write.csv(data, "<file_1.csv>")

#############

conditions = factor(sel_samples$MMR_IHC)

pvalues_all <- sapply(1:nrow(normcounts),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(normcounts[i,])),conditions)
  wtest=wilcox_test(gene~conditions, data, distribution = "exact")
  p = pvalue(wtest)
  return(p)
})

fdr=p.adjust(pvalues_all,method = "fdr")

d_mmr = normcounts[,sel_samples$MMR_IHC == 'Deficient']
p_mmr = normcounts[,sel_samples$MMR_IHC == 'Proficient']
log2fc = log2(rowMeans(d_mmr) / rowMeans(p_mmr))
de_genes = data.frame(log2FoldChange = log2fc, pValues = pvalues_all, adjusted_pValues = fdr)
rownames(de_genes) = rownames(normcounts)

library("biomaRt")
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')

annot <- getBM(
  attributes = c(
    'hgnc_symbol',
    'external_gene_name',
    'ensembl_gene_id_version',
    'gene_biotype',
    'description'
  ),
  filters = 'ensembl_gene_id_version',
  values = rownames(normcounts),
  mart = ensembl)


annot <- merge(
  x = de_genes,
  y =  annot,
  by.y = 'ensembl_gene_id_version',
  all.x = T,
  by.x = "row.names")



annot$diffexpressed <- "NO"
# if log2Foldchange > 1 and ajusted_pvalue < 0.05, set as "UP" 
annot$diffexpressed[annot$log2FoldChange > 1 & annot$adjusted_pValues < 0.05] <- "UP"
# if log2Foldchange < -1 and ajusted_pvalue < 0.05, set as "DOWN"
annot$diffexpressed[annot$log2FoldChange < -1 & annot$adjusted_pValues < 0.05] <- "DOWN"

result = annot[annot$diffexpressed != "NO",]
result = result[order(result$log2FoldChange),]
write.csv(result, "<file_2.csv>")


