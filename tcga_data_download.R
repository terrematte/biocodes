#'-----------------------------------------------------------------
# RTCGAToolbox ----
#'-----------------------------------------------------------------
packages_bioconductor = c("SummarizedExperiment", "RTCGAToolbox")
if (!require(x, character.only = TRUE)) {
  BiocManager::install(x, dependencies = TRUE)
  library(x, character.only = TRUE)
}

# Show all datasets
getFirehoseDatasets()

# Show the last three updates
getFirehoseRunningDates(last = 3)

# Get mutation data and clinical data
kircData <- getFirehoseData(dataset="KIRC", runDate="20160128", miRNASeqGene = T, RNASeqGene = T,
                            forceDownload=TRUE, clinical=TRUE, Mutation=F)

# Select RNASeq data
kirc.rnaseq <- biocExtract(kircData, "RNASeqGene")

# Select RNASeq data
kirc.exp <- assay(kirc.rnaseq)
kirc.exp <- as.data.frame(kirc.exp)

kirc.exp["AAAS", c("TCGA-A3-3306-01A-01R-0864-07", "TCGA-A3-3307-01A-01R-0864-07", "TCGA-A3-3308-01A-02R-1325-07")]

#'-----------------------------------------------------------------
# Xenabrowser ----
#'-----------------------------------------------------------------

# From https://xenabrowser.net/datapages/?dataset=TCGA-KIRC.htseq_counts.tsv&host=https%3A%2F%2Fgdc.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

url <- "https://gdc-hub.s3.us-east-1.amazonaws.com/latest/TCGA-KIRC.htseq_counts.tsv.gz"
destfile <- "kirc_counts.tsv.gz"
download.file(url, destfile)
library(tidyverse)
kirc.counts <- read_tsv(gzfile("kirc_counts.tsv.gz"))

# Select only the Ensembl, an remove the version
#row.names(kirc.counts) <- sapply(strsplit(kirc.counts$Ensembl_ID, "\\."), `[`, 1)
row.names(kirc.counts) <- sub("\\..*", "", kirc.counts$Ensembl_ID)
kirc.counts$Ensembl_ID <- NULL
kirc.counts <- as.data.frame(kirc.counts)

# Revert the normalization log2(count + 1)
kirc.counts <- 2^(kirc.counts)-1

rm(url, destfile)

#'-----------------------------------------------------------------
# Firehose ----
#'-----------------------------------------------------------------

# From Firehose
url <- " http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/KIRC/20160128/gdac.broadinstitute.org_KIRC.Merge_rnaseq__illuminahiseq_rnaseq__unc_edu__Level_3__gene_expression__data.Level_3.2016012800.0.0.tar.gz"
destfile <- "kirc_counts.tar.gz"
download.file(url, destfile)
library(tidyverse)

KIRC.rnaseq <- read.delim(untar("kirc_counts.tar.gz"))

KIRC.rnaseq.count <- KIRC.rnaseq %>% 
  column_to_rownames("Hybridization.REF")  %>% 
  dplyr::select(ends_with(".07"))


#'-----------------------------------------------------------------
# TCGABiolinks ----
#'-----------------------------------------------------------------
# Installation from GitHub
# BiocManager::install("BioinformaticsFMRP/TCGAbiolinks")
#' https://www.bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
packages_bioconductor = c("TCGAbiolinks", "SummarizedExperiment")
package.check <- lapply(packages_bioconductor, FUN = function(x) {
  if (!require(x, character.only = TRUE)) {
    BiocManager::install(x, dependencies = TRUE)
    library(x, character.only = TRUE)
  }
})

query <- GDCquery(project = "TCGA-KIRC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "HTSeq - Counts",
                  sample.type = c("Primary solid Tumor", "Solid Tissue Normal"))
biomartCacheClear() 
GDCdownload(query, files.per.chunk = 100)
data <- GDCprepare(query)