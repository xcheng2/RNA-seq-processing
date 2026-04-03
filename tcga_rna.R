rm(list = ls())
setwd("<dir_name>")

library(TCGAbiolinks)

query <- GDCquery(
  project = c("TCGA-COAD", "TCGA-READ"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
#  workflow.type = "STAR - Counts",
#  sample.type = "Primary Tumor"
)


coad_samples = getResults(query)

GDCdownload(query, directory = "./TCGA_COAD_READ")

coad_data <- GDCprepare(query, directory = "./TCGA_COAD_READ")

#######
raw_count = coad_data@assays@data$unstranded

raw_count = as.data.frame(raw_count, row.names = rownames(coad_data))

colnames(raw_count) = coad_samples$cases

tpm = coad_data@assays@data$tpm_unstrand

tpm = as.data.frame(tpm, row.names = rownames(coad_data))

colnames(tpm) = coad_samples$cases


  

######
# remove some duplicate
table(coad_samples$sample_type)

#Metastatic       Primary Tumor     Recurrent Tumor Solid Tissue Normal 
#1                 647                   2                  51 
normal = coad_samples[coad_samples$sample_type == 'Solid Tissue Normal',]
t_n = table(normal$cases.submitter_id)
max(t_n)

# no duplicate for normal tissue

sel = coad_samples[coad_samples$sample_type == "Primary Tumor" | coad_samples$sample_type == "Solid Tissue Normal", ]
tumor = sel[sel$sample_type == "Primary Tumor", ]

table(substr(tumor$cases, 20, 20))
# Follow Broad Institute's rule on selecting from multiple
#https://broadinstitute.atlassian.net/wiki/spaces/GDAC/pages/844334036/FAQ#%5BinlineExtension%5DQ%3A-What-do-you-do-when-multiple-aliquot-barcodes-exist-for-a-given-sample%2Fportion%2Fanalyte-combination%3F
# all are R (H > R > T), so prefer the highest lexicographical sort value 
tumor = tumor[order(tumor$cases, decreasing = TRUE), ]
rownames(tumor) = tumor$cases

temp = ""
mylst = c()
for (i in 1:nrow(tumor)) {
    if (tumor$cases.submitter_id[i] == temp) mylst = c(mylst, list(tumor$cases[i]))
    temp = tumor$cases.submitter_id[i]
}


final_sample = sel[!sel$cases %in% mylst, ]
write_csv(final_sample, "TCGA_unique_sample.csv")
save.image(file = "tcga_crc.RData")

################

load("tcga_crc.RData")
new_tumor = final_sample[final_sample$sample_type == "Primary Tumor", ]

library(edgeR)

sel_count = raw_count[,new_tumor$cases]

colnames(sel_count) = substr(colnames(sel_count), 1, 12)
sel_count$summ = rowSums(sel_count[, 1:624])
sel_count = sel_count[sel_count$summ > 0,]

rownames(sel_count) = substr(rownames(sel_count), 1, 15)
y = DGEList(counts = sel_count[, 1:624], genes = rownames(sel_count))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method="TMM")
normcounts<-cpm(y, normalized.lib.sizes = T)


CXCL10 = log2(normcounts['ENSG00000169245',]+1)
CCL5 = log2(normcounts['ENSG00000271503',]+1)
CXCL9 = log2(normcounts['ENSG00000138755',] +1)
CXCL11 = log2(normcounts['ENSG00000169248',] +1)
CTNNB1 = log2(normcounts['ENSG00000168036',] + 1)
GREM1 = log2(normcounts['ENSG00000166923',] + 1)

data = data.frame(CXCL10_log2normCounts = CXCL10, CCL5_log2normCounts = CCL5, CXCL9_log2normCounts = CXCL9, CXCL11_log2normCounts = CXCL11, CTNNB1_log2normCounts = CTNNB1, GREM1_log2normCounts = GREM1)

write.csv(merged, "<filename.csv>")
