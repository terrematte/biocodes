
gene_missing <- c("C10orf4", "C17orf79", "NCRNA00116", "ODZ3", "Hsa-mir-374a", "ERVFRDE1", "FAM123A",
                  "VEGFR1", "VEGFR2", "VEGFR3", "HIF2A", "PD-L1", "IFN-G" , "PD1", "PD-L2", "TIM3" , "CNN1G",
                  "EXOC3L2", "FYB" ,"ADRBK1", "FAM178A", "BAT2")

#' -----------------------------------------------------
# org.Hs.eg.db ----
#' -----------------------------------------------------

packages_bioconductor = c("AnnotationDbi", "org.Hs.eg.db")
if (!require(x, character.only = TRUE)) {
  BiocManager::install(x, dependencies = TRUE)
  library(x, character.only = TRUE)
}


# Check keytypes available
columns(org.Hs.eg.db)

#To map one key to a column
#gene_miss_df <- mapIds(org.Hs.eg.db, keys=gene_missing, column = "SYMBOL", keytype=c("ENSEMBL","SYMBOL", "GENENAME", "ONTOLOGY"))
cols <- c("ENSEMBL","SYMBOL", "GENENAME")
gene_miss_df <- select(org.Hs.eg.db, keys=gene_missing, columns=cols, keytype="SYMBOL")
gene_miss_df

#' -----------------------------------------------------
# mygene ----
#' -----------------------------------------------------

packages_bioconductor = c("mygene")
if (!require(x, character.only = TRUE)) {
  BiocManager::install(x, dependencies = TRUE)
  library(x, character.only = TRUE)
}

gene_miss_df <- queryMany(gene_missing, scopes="symbol", fields = c("symbol", "name", "taxid"), species="human")
gene_miss_df

#' -----------------------------------------------------
# biomaRt ----
#' -----------------------------------------------------

library("biomaRt")
library(ensembldb)
library(AnnotationHub)
ensembl <- useMart("ensembl")

ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)

gene_miss_df <- biomaRt::getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'description'), 
                               filters = 'hgnc_symbol', 
                               values = gene_missing, 
                               mart = ensembl,  
                               useCache = FALSE)
#' -----------------------------------------------------
## Alterative way to get gene Annotations by biomart package ----
#' -----------------------------------------------------

ennsembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

ensembl_genes <- row.names(exp)

ensemblAnnot <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id", "description"),
                      filters=c("hgnc_symbol"),
                      values=ensembl_genes, mart=ensembl)


genes_miss <- setdiff(ensembl_genes, unique(ensemblAnnot$hgnc_symbol))

while(length(genes_miss) > 10){
  
  for(i in seq(from=1, to=length(genes_miss), by=100)){
    try(ensemblAnnot_1 <- getBM(attributes= c("hgnc_symbol", "ensembl_gene_id", "description"),
                                filters=c("hgnc_symbol"),
                                values=genes_miss[i:(i+99)], mart=ensembl), silent = T)
    ensemblAnnot <- rbind(ensemblAnnot, ensemblAnnot_1)
    ensemblAnnot <- ensemblAnnot[!duplicated(ensemblAnnot),]
    genes_miss <- setdiff(ensembl_genes, unique(ensemblAnnot$hgnc_symbol))
    print(i)
    Sys.sleep(1)
  }
  genes_miss <- setdiff(ensembl_genes, unique(ensemblAnnot$hgnc_symbol))
}

colnames(ensemblAnnot) <- c("GeneID", "Ensembl", "Description")

rm(ensembl,ensembl_genes)

#' -----------------------------------------------------
# annotables ----
#' -----------------------------------------------------

#' If you haven't already installed devtools...
# install.packages("devtools")
#'
#' Use devtools to install the package
#devtools::install_github("stephenturner/annotables")

library(annotables)

ensemblAnnot <- grch38 %>%
  filter(grch38$symbol %in%  gene_missing)  %>%
  dplyr::select(ensgene, symbol, description)

#' -----------------------------------------------------
# GeneBook ----
#' -----------------------------------------------------

#install.packages("GeneBook")
library("GeneBook")

##  Multiple Gene ID Convert
gene_missing_m <- as.matrix(gene_missing)
mat_id_convert = c()
for(i in 1:nrow(gene_missing_m)){
  out <- GeneCard_ID_Convert(gene_missing_m[i])
  mat_id_convert=rbind(mat_id_convert,out)
}

gene_missing_results <- cbind(gene_missing_m, mat_id_convert)
colnames(gene_missing_results) <- c("previous_ID","Symbol","Label")
gene_missing_results <- data.frame(gene_missing_results)
