library(TCGAbiolinks)
library(readxl)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(MultiAssayExperiment)
library(ELMER.data)
library(ELMER)
library(dplyr)
library(plyr)
library(png)
# devtools::install_version("dbplyr", version = "2.3.4")
plyr

dir.data <- "data"
dir.plot <- "plots"
dir.tcga <- file.path(dir.data, "TCGA")
dir.output <- "analysis_results/TCGA/ELMER_analysis_hg19"
for (p in grep("dir", ls(), value = T)) dir.create(get(p), recursive = TRUE, showWarnings = FALSE)

metadata <- readxl::read_xlsx(
  path = file.path(dir.data, "Table_S2_TNBCsubtype clinical information and signatures.xlsx"), 
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS")

load(file.path(dir.tcga, "TCGA_TNBC_192samples_RNASeq.Rdata"))
rna <- dataFilt

rownames(rna) <- 
  rowRanges(BRCA.exp)$gene_id[
    match(rownames(rna),
          (adply(rowRanges(BRCA.exp)$gene_id, .margin = 1, .fun = function(s){
            y = unlist(strsplit(s, "[.]"))[1]
            return(y)
            }))[,2]
          )]

rna <- rna[!is.na(rownames(rna)),]



met <- get(load(file.path(dir.tcga, "TCGA_TNBC_192samples_Methylation.Rdata")))
genome <- "hg19"
met.platform <- "450K"
distal.probes <- get.feature.probe(genome = genome, met.platform = met.platform)

mae.distal <- createMAE(exp = rna, met = met)

