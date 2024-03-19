# load data
load("data/TCGA/TCGA_TNBC_192samples_RPPA.Rdata")

exp <- BRCA.RPPA
rppa.samples <- substr(colnames(BRCA.RPPA), 1, 12)[-1] %>% unique
rownames(exp) <- exp$peptide_target
exp <- exp[!is.na(rowSums(exp[, -(1:5)])), -(1:5)]
rownames(exp)

exp[]

exp[, -(1:5)]
substr(colnames(BRCA.RPPA), 1, 12)[-1] %>% unique
rownames(BRCA.RPPA)
dim(exp)
sum(BRCA.RPPA$catalog_number == "CTCF")
