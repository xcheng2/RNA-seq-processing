rm(list = ls())

setwd("<dir_name>")

library(edgeR)

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

conditions = factor(sel_samples$Race)

design <- model.matrix(~conditions)
y <- estimateDisp(y, design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit, coef=2)
de_genes = as.data.frame(topTags(qlf, n=nrow(y))$table)
rownames(de_genes) = de_genes$Geneid

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
  values = rownames(de_genes),
  mart = ensembl)


annot <- merge(
  x = de_genes,
  y =  annot,
  by.y = 'ensembl_gene_id_version',
  all.x = T,
  by.x = "row.names")

annot$logFC = -annot$logFC

annot$diffexpressed <- "NO"
# if log2Foldchange > 1 and ajusted_pvalue < 0.05, set as "UP" 
annot$diffexpressed[annot$logFC > 1 & annot$FDR < 0.05] <- "UP"
# if log2Foldchange < -1 and ajusted_pvalue < 0.05, set as "DOWN"
annot$diffexpressed[annot$logFC < -1 & annot$FDR < 0.05] <- "DOWN"

result = annot[annot$diffexpressed != "NO",]

result = result[order(result$logFC),]
write.csv(result, "<out_filename.csv>")
