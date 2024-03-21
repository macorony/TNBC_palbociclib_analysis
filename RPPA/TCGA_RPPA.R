# load data
library(limma)
library(edgeR)
library(plyr)
library(sva)
library(dplyr)
options(scipen = 500)

dir.data <- "data"
dir.plots <- "plots"
dir.tcga <- file.path(dir.data, "TCGA")
dir.output <- "analysis_results/TCGA/rppa_results"
for (p in grep("dir", ls(), value = T)) dir.create(get(p), recursive = TRUE, showWarnings = FALSE)

load("data/TCGA/TCGA_TNBC_192samples_RPPA.Rdata")

rppa.samples <- substr(colnames(BRCA.RPPA)[-(1:5)], 1, 12) %>% unique

BRCA.RPPA.avg <- plyr::adply(rppa.samples, .margins = 1, .fun = function(s){
  rowMeans(BRCA.RPPA[, grep(s, colnames(BRCA.RPPA)), drop = FALSE], na.rm = TRUE)
  }, .id = NULL) %>% t
rownames(BRCA.RPPA.avg) <- BRCA.RPPA$peptide_target
colnames(BRCA.RPPA.avg) <- substr(colnames(BRCA.RPPA),1,15)[-(1:5)] %>% unique

exp <- BRCA.RPPA.avg
exp <- exp[!is.na(rowSums(exp)), ]
colnames(exp) <- substr(colnames(exp), 1, 12)

metadata <- readxl::read_xlsx(path = file.path(dir.data, "Table_S2_TNBCsubtype clinical information and signatures.xlsx"), 
                              sheet = "A-TCGA_TNBC_subtype") %>% 
  dplyr::filter(subtype != "UNS") %>% 
  dplyr::filter(!is.na(Stage))

exp <- exp[, colnames(exp) %in% metadata$patient]
metadata <- metadata[match(colnames(exp), metadata$patient), ]

DM_limma <- function(
    data, 
    sample_info,
    grp_by,
    mult_factor, 
    output_dir = NULL, 
    prefix = NULL, 
    output_name = NULL){
  Temp_id <- na.omit(match(colnames(data), sample_info$patient))
  sample_info <- sample_info[Temp_id, ]
  Temp_id <- na.omit(match(sample_info$patient, colnames(data)))
  data <- data[, Temp_id]
  cat(paste0("~ 0 + ", grp_by, " + Stage + Age"), "\n")
  sample_info$Stage <- factor(sample_info$Stage)
  mod <- model.matrix(formula(paste0("~", mult_factor)), data = sample_info)
  mod0 <- model.matrix(~Age+Stage, data = sample_info)
  
  data <- data %>% na.omit %>% as.matrix()
  
  n.sv <- num.sv(data, mod, method = "be", B =20, seed = set.seed(5000))
  n.sv
  
  svobj <- sva(data, mod, mod0, n.sv = n.sv)
  
  df <- svobj$sv %>% as.data.frame
  colnames(df) <- paste0("sv", 1:n.sv)
  

}
metadata$Stage = factor(metadata$Stage)
model.matrix(formula(paste0("~", "BL1+Stage+Age")), data = metadata)
model.matrix(~Age+Stage, data= metadata)
xnam <- paste0("x", 1:25)
as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
paste("y ~ ", paste(xnam, collapse = "+"))


