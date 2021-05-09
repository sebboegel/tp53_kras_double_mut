library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)
library(dplyr)
library(biomaRt)

#retrieve mutation data (hg38)
#https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/mutation.html)
maf <- GDCquery_Maf("LUAD", pipelines = "muse")
maf_tp53=filter(maf,Hugo_Symbol == "TP53")
maf_kras=filter(maf,Hugo_Symbol == "KRAS")
tp53_kras_mut=intersect(maf_tp53$Tumor_Sample_Barcode,maf_kras$Tumor_Sample_Barcode)#barcodes


tp53=as.data.frame(maf_tp53)
kras=as.data.frame(maf_kras)
tp53_kras=rbind(tp53,kras)
tp53_kras_overlap=tp53_kras[tp53_kras$Tumor_Sample_Barcode %in% tp53_kras_mut,]
write.csv2(tp53_kras_overlap,"tp53_kras_mut_info_overlap.csv")
tp53_kras_mut=substr(tp53_kras_mut,0,16)#patient_ids

#patients not having neither a tp53 nor a kras mutation
maf_tp53_not=filter(maf,Hugo_Symbol != "TP53")
maf_kras_not=filter(maf,Hugo_Symbol != "KRAS")
tp53_kras_mut_not=intersect(maf_tp53_not$Tumor_Sample_Barcode,maf_kras_not$Tumor_Sample_Barcode)#barcodes

tp53_not=as.data.frame(maf_tp53_not)
kras_not=as.data.frame(maf_kras_not)
tp53_kras_not=rbind(tp53_not,kras_not)
tp53_kras_not_overlap=tp53_kras_not[tp53_kras_not$Tumor_Sample_Barcode %in% tp53_kras_mut_not,]
write.csv2(tp53_kras_not_overlap,"not_tp53_kras_mut_info_overlap.csv")

tp53_kras_mut_not=substr(tp53_kras_mut_not,0,16)#patient_ids
tp53_kras_mut=substr(tp53_kras_mut,0,16)#patient_ids
#retrive all available RNA-Seq barcodes for LUAD
## Gene expression aligned against hg38
#https://bioconductor.org/packages/devel/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
query <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification", workflow.type = "HTSeq - FPKM-UQ")
GDCdownload(query)
data <- GDCprepare(query)
barcodes=as.data.frame(assay(data))
barcodes=names(barcodes)

#double mutated
list_intersect_barcodes=sapply(tp53_kras_mut, function(x) barcodes[grepl(x, barcodes)])
list_intersect_barcodes=unlist(list_intersect_barcodes,use.names = F)

#retrieve all gene expression values for intersect barcodes
query <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "HTSeq - FPKM-UQ", barcode = list_intersect_barcodes)
GDCdownload(query)
data <- GDCprepare(query)
my_gene_exp=as.data.frame(assay(data))
my_gene_exp <- cbind(rownames(my_gene_exp), data.frame(my_gene_exp, row.names=NULL))
gene_exp_mutated=rbind(my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000136352"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000149948"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000149021"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000168484"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000181449"),])
write.csv2(gene_exp_mutated,"gene_exp_double_mutated.csv")

#not mutated
list_intersect_barcodes=sapply(tp53_kras_mut_not, function(x) barcodes[grepl(x, barcodes)])
list_intersect_barcodes=unlist(list_intersect_barcodes,use.names = F)

#retrieve all gene expression values for intersect barcodes
query <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling",data.type = "Gene Expression Quantification",workflow.type = "HTSeq - FPKM-UQ", barcode = list_intersect_barcodes)
GDCdownload(query)
data <- GDCprepare(query)
my_gene_exp=as.data.frame(assay(data))
my_gene_exp <- cbind(rownames(my_gene_exp), data.frame(my_gene_exp, row.names=NULL))
gene_exp_mutated=rbind(my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000136352"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000149948"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000149021"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000168484"),],my_gene_exp[which(my_gene_exp$`rownames(my_gene_exp)`=="ENSG00000181449"),])
write.csv2(gene_exp_mutated,"gene_exp_not_mutated.csv")
