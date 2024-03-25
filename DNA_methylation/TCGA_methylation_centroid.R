library(Hmisc)
library(qvalue)
library(limma)
library(edgeR)
library(hexbin)
library(sva)
library(dplyr)
require(SummarizedExperiment)
options(scipen = 500)

dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data, "TCGA")
dir.output <- "analysis_results/TCGA/dnam_results"
for (p in grep("dir", ls(), value = T)) dir.create(get(p), recursive = T, showWarnings = F)

load(file.path(dir.tcga, "TCGA_TNBC_192samples_Methylation.Rdata"))
dnam <- assay(BRCA.met)
colnames(dnam) <- substr(colnames(dnam), 1, 12)

BMIQ.nonrs <- dnam[!grepl("rs.", rownames(dnam)),] %>% na.omit()

BMIQ <- log2(BMIQ.nonrs/(1 - BMIQ.nonrs))

metadata <- readxl::read_xlsx(
  path = file.path(dir.data, "Table_S2_TNBCsubtype clinical information and signatures.xlsx"),
  sheet = "A-TCGA_TNBC_subtype"
) %>% dplyr::filter(subtype != "UNS") %>% dplyr::filter(!is.na(Stage))

BMIQ <- BMIQ[, colnames(BMIQ) %in% metadata$patient]
metadata <- metadata[match(colnames(BMIQ), metadata$patient),]

DM_limma <- function(
    data, 
    sample_info,
    grp_by,
    mult_factor, 
    output_dir = NULL, 
    prefix = NULL){
  
  stopifnot(all(colnames(data) == sample_info$patient))
  cat(paste0("~ 0 +", grp_by, " + Stage + Age"), "\n")
  
  # preparing the data
  sample_info$Stage = factor(sample_info$Stage)
  mod <- model.matrix(formula(paste0("~ 0 +", mult_factor)), data = sample_info)
  mod0 <- model.matrix(~ Age + Stage, data = sample_info)
  
  data <- data %>% na.omit %>% as.matrix
  
  n.sv <- num.sv(data, mod, method = "be", B = 20, seed = set.seed(5000))
  n.sv
  
  svobj <- sva(data, mod, mod0, n.sv = n.sv)
  
  df <- svobj$sv %>% as.data.frame
  colnames(df) <- paste0("sv", 1:n.sv)
  modSv = cbind(mod, df)
  mod0Sv = cbind(mod0, df)
  
  data.vlmfit = limFit(data, modSv)
  Fit.eb <- eBayes(data.vlmfit)
}

