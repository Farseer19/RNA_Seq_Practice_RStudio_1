install.packages("BiocManager")
BiocManager::install("recount3")
BiocManager::install("VariantExperiment")
BiocManager::install("DESeq2")
BiocManager::install("apeglm")

library("VariantExperiment")
library("recount3")
library("DESeq2")
library("stringi")
library("apeglm")




rse_gene_SRP126485 <- recount3::create_rse_manual(
  project = "SRP126485",
  project_home = "data_sources/sra",
  organism = "human",
  annotation = "gencode_v26",
  type = "gene"
)

rse_gene_SRP126485

#Creating a list of all sample titles to create columns for DESeq2 design
colN <- rse_gene_SRP126485@colData@listData@raw_counts[28]

#Vectors that will identify patient and treatment type
patient <- c()
treatment <- c()


#Loop through each sample title to extract patient and treatment data and add to vector
for (index in colN[[1]]) {
  if (grepl("Placebo",index)){
    treatment <- c(treatment, "Placebo")
  }else{
    treatment <- c(treatment, "Metformin")
  }
  patient <- c(patient, stri_extract_first_regex(index, "[0-9]+"))
}


patient
treatment
#Add new vectors to colData
rse_gene_SRP126485@colData@listData[["patient"]] <- as.integer(patient)
rse_gene_SRP126485@colData@listData[["treatment"]] <- treatment

#Create DESeq Data Set from RSE Object
rse_gene <- DESeqDataSet(rse_gene_SRP126485, design = ~ patient + treatment)
rse_gene <- DESeq(rse_gene)

#Check to see non-zero counts in assay
colA <- rse_gene_SRP126485@assays@data@listData[[1]]
colA <- colA[!(colA %in% 0)]
colA

#Pre-filter data before analysis - Ensures that there will be a count of at least 10 for each row and minimum number of samples will be set to 58 as that is the number of treated samples
smallestGroupSize <- 58
keep <- rowSums(counts(rse_gene) >= 10) >= smallestGroupSize
rse_gene <- rse_gene[keep,]



#Generate Results from pre-filtered data
res <- results(rse_gene)
res
summary(res)

#Log Fold Change Shrinkage
resLFC <- lfcShrink(rse_gene, coef = "treatment_Placebo_vs_Metformin", type="apeglm")
resLFC
summary(resLFC)

#Find genes that are relevant through filtering for a lower adjusted p-value
res_subset <- subset(resLFC, padj < 0.2)
res_subset

#Plotting both the shrunken results and the noisy
plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))






print(citation("recount3")[2], bibtex = TRUE) 

