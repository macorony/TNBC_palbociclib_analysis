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

mae.distal <- createMAE(exp = rna, met = met, save = FALSE, linearize.exp = TRUE, 
                        filter.probes = distal.probes, met.platform = met.platform, 
                        genome = genome, TCGA = TRUE)
mae.distal$TNBC_subtype <- gsub("[0-9]-", "", 
                                metadata$subtype[match(mae.distal$patient, metadata$patient)])

for (i in unique(mae.distal$TNBC_subtype)) {
  colData(mae.distal)[[paste0(i, "_vs_others")]] <- "Others"
  colData(mae.distal)[[paste0(i, "_vs_others")]][mae.distal$TNBC_subtype == i] <- i
}

save(mae.distal, file = "data/mae_TNBC_with_groups_distal_probes_hg19.rda")

mae.distal = mae.distal[, !is.na(mae.distal$TNBC_subtype)]

cores <- 1
mode <- "supervised"

for(tnbc.subtype in unique(mae.distal$TNBC_subtype)) {
  for(direction in c(c("hypo", "hyper"))) {
    group.col <- paste0(tnbc.subtype, "_vs_others")
    group1 <- tnbc.subtype
    group2 <- "Others"
    message(group.col, ":", direction, " ", group1, " vs ", group2)
    
    dir.out <- paste0(
      dir.output, "/", mode, "/", group.col, "-", 
      gsub("[[:punct:]]| ", ".", group1), 
      "_vs_", 
      gsub("[[:punct:]]| ", ".", group2),
      "/", direction
    )
    sig.diff <- get.diff.meth(
      data = mae.distal,
      group.col = group.col,
      group1 = group1,
      group2 = group2, 
      mode = mode,
      sig.dif = 0.2,
      diff.dir = direction, 
      cores = cores, 
      dir.out = dir.out,
      pvalue = 0.01
    )
    nearGenes <- GetNearGenes(
      data = mae.distal, 
      probes = sig.diff$probe, 
      numFlankingGenes = 20
    )
    pair <- get.pair(
      data = mae.distal, 
      group.col = group.col,
      group1 = group1,
      group2 = group2,
      nearGenes = nearGenes, 
      mode = mode, 
      permu.dir = paste0(dir.out, "/permu"),
      permu.size = 10000, 
      raw.pvalue = 0.05,
      Pe = 0.05,
      filter.probes = T,
      filter.percentage = 0.05,
      filter.portion = 0.3,
      dir.out = dir.out,
      diff.dir = direction,
      cores = cores,
      label = direction
    )
    tryCatch({
      heatmapPairs(
        data = mae.distal,
        group.col = group.col,
        group1 = group1,
        group2 = group2,
        annotation.col = c("TNBC_subtype"),
        pairs = pair,
        filename = paste0(dir.out, "/heatmap_pairs.pdf")
    )
      }, error = function(e) print(e))
    
    if(length(unique(pair$Probe)) < 10) next
    
    enriched.motif <- get.enriched.motif(
      data = mae.distal, 
      probes = unique(pair$Probe), 
      dir.out = dir.out, 
      label = direction,
      min.incidence = 10, 
      lower.OR = 1.1
    )
    
    TF <- get.TFs(
      data = mae.distal, 
      group.col = group.col,
      group1 = group1,
      group2 = group2,
      mode = mode,
      diff.dir = direction,
      enriched.motif = enriched.motif,
      dir.out = dir.out,
      cores = cores,
      label = direction,
      save.plots = F
    )
  
  }
}

