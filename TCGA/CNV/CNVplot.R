library(dplyr)
library(TCGAbiolinks)
library(GenomicRanges)
library(ggpubr)


# How to create folders
dir.data <- "data/GISTIC/"
dir.plot <- "palbociclib_plots/GISTIC"
dir.analysis <- "palbociclib_analysis/TCGA/CNV"
dir.output <- "palbociclib_analysis_results/GISTIC"

for (p in grep("dir", ls(), value = T)) dir.create(get(p), recursive = T, showWarnings = FALSE)

# Input the necessary function
source("analysis/TCGA/CNV/VisualizeCNplot.R")
# Gene information
genes <- TCGAbiolinks::get.GRCh.bioMart("hg19") %>%
  dplyr::select(external_gene_name, chromosome_name, start_position, end_position) %>%
  dplyr::rename(GeneSymbol = external_gene_name, Chr = chromosome_name, Start = start_position, End = end_position) %>% 
  unique()

GeneToVisualize <- c("RSL24D1", "MAK16", "ARHGAP19", "CNOT6L", "CTCF", "SFPQ", "UTP23", "HNRNPA1", "CYCS", "RBMX", "RFWD3", "SENP1",
                     "TMA16", "NRF1", "DIMT1", "LIN54", "RFX7", "DPH5", "SKIV2L2", "SNRNP40", "C7orf26", "TOMM22", "IRF2BP2", "RIOK2",
                     "MARS2", "SERBP1", "DTL", "INTS7", "UBTF", "RBM8A", "WDR43", "CREBBP", "HNRNPF", "PCBP2", "RFXAP", "NASP", "PM20D2",
                     "WDR82", "NOP14", "POLR1B", "TIMELESS", "NCAPD3", "BRD2")
GeneToVisualize <- sort(unique(GeneToVisualize))

FolderList <- list.dirs(dir.data)
# don't write function
curSubtype <- "TNBCall"
folderCur <- FolderList[grep(curSubtype, FolderList)]
FolderListCur <- list.files(folderCur)

fileCurScore <- FolderListCur[grep("scores.gistic", FolderListCur)]
fileCurAmp <- FolderListCur[grep("amp_genes.conf_99", FolderListCur)]
fileCurDel <- FolderListCur[grep("del_genes.conf_99", FolderListCur)]

tabAnno_conf_Amp <- readr::read_delim(paste0(folderCur, "/", fileCurAmp), "\t", escape_double = FALSE, trim_ws = TRUE)
tabAnno_conf_Del <- readr::read_delim(paste0(folderCur, "/", fileCurAmp), "\t", escape_double = FALSE, trim_ws = TRUE)
scores_filt_call <- readr::read_delim(paste0(folderCur, "/",fileCurScore), "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  as.data.frame %>% dplyr::mutate(Aberration.Kind = 0) %>% 
  dplyr::mutate(Aberration.Kind = replace(Type, Type == "Amp", 1)) %>%
  dplyr::mutate(Aberration.Kind = replace(Type, Type == "Del", -1)) %>%
  dplyr::rename("Region Start [bp]" = Start, "Region End [bp]" = End, "score" = "G-score") %>%
  dplyr::mutate(Ratio = score/`-log10(q-value)`)

FDR.thresh <- 1 / 10 ^ 0.25

scores_filt_call_Amp_full <- scores_filt_call %>% dplyr::filter(Type == "Amp")
scores_filt_call_Del_full <- scores_filt_call %>% dplyr::filter(Type == "Del")
scores_filt_FDR_Amp <- scores_filt_call_Amp_full %>% dplyr::filter(`-log10(q-value)` > -log10(FDR.thresh))
scores_filt_FDR_Del <- scores_filt_call_Del_full %>% dplyr::filter(`-log10(q-value)` > -log10(FDR.thresh))

genes_GR <- makeGRangesFromDataFrame(genes, keep.extra.columns = TRUE)

# Adding gene names to the gistic result
df_GR <- scores_filt_call_Amp_full
colnames(df_GR)[colnames(df_GR) == "Region Start [bp]"] <- "start"
colnames(df_GR)[colnames(df_GR) == "Region End [bp]"] <- "end"
df_GR_new <- makeGRangesFromDataFrame(df_GR, keep.extra.columns = TRUE)
hits_Amp <- findOverlaps(genes_GR, df_GR_new, type = 'any')
df_GR_ann_Amp <- cbind(df_GR[subjectHits(hits_Amp),], genes[queryHits(hits_Amp),])
df_GR_ann_Amp <- df_GR_ann_Amp[order(df_GR_ann_Amp$score, decreasing = TRUE),]
df_GR_ann_Amp <- df_GR_ann_Amp[!duplicated(df_GR_ann_Amp$GeneSymbol),]

GenesinWidePeak_Amp <- NULL
if (ncol(tabAnno_conf_Amp) > 1) {
  # search for amplified genes
  for (curColumn in grep("[[:alnum:]]*p|[[:alnum:]]q", colnames(tabAnno_conf_Amp), value = T)) {
    curGeneWidePeak <- intersect(as.matrix(tabAnno_conf_Amp[, curColumn]), df_GR_ann_Amp$GeneSymbol)
    print(sort(curGeneWidePeak))
    GenesinWidePeak_Amp <- c(curGeneWidePeak, GenesinWidePeak_Amp)
  }
  df_GR_ann_Amp_GeneWidePeak <- df_GR_ann_Amp %>% 
    dplyr::filter(GeneSymbol %in% GenesinWidePeak_Amp)
} else if (ncol(tabAnno_conf_Amp) == 1) {
  df_GR_ann_Amp_GeneWidePeak <- tabAnno_conf_Amp
}

df_GR <- scores_filt_call_Del_full
colnames(df_GR)[colnames(df_GR) == "Region Start [bp]"] <- "start"
colnames(df_GR)[colnames(df_GR) == "Region End [bp]"] <- "end"
df_GR$Chromosome <- gsub("chr0", "", df_GR$Chromosome)
df_GR$Chromosome <- gsub("chr","", df_GR$Chromosome)
df_GR$Chromosome <- as.numeric(df_GR$Chromosome)
df_GR <- df_GR[df_GR$Type %in% "Del", ]
df_GR_new <- makeGRangesFromDataFrame(df_GR, keep.extra.columns = TRUE)
hits_Del <- findOverlaps(genes_GR, df_GR_new, type = "any")
df_GR_ann_Del <- cbind(df_GR[subjectHits(hits_Del), ], genes[queryHits(hits_Del), ])
df_GR_ann_Del <- df_GR_ann_Del[order(df_GR_ann_Del$score, decreasing = TRUE),]
df_GR_ann_Del <- df_GR_ann_Del[!duplicated(df_GR_ann_Del$GeneSymbol), ]

GenesinWidePeak_Del <- NULL
if (ncol(tabAnno_conf_Del) >1) {
  for (curColumn in grep("[[:alnum:]]*p|[[:alnum:]]q",colnames(tabAnno_conf_Del),value = T)) {
    curGeneWidePeak <- intersect(as.matrix(tabAnno_conf_Del[, curColumn]), df_GR_ann_Del$GeneSymbol)
    print(sort(curGeneWidePeak))
    GenesinWidePeak_Del <- c(curGeneWidePeak, GenesinWidePeak_Del)
  }
  
  df_GR_ann_Del_GeneWidePeak <- df_GR_ann_Del %>%
    dplyr::filter(GeneSymbol %in% GenesinWidePeak_Del)
} else if (ncol(tabAnno_conf_Del) == 1) {
  df_GR_ann_Del_GeneWidePeak <- tabAnno_conf_Del
}

transformationFormulaAxis2 <- rbind(scores_filt_FDR_Amp, scores_filt_FDR_Del)
transformationFormulaAxis2_thresh <- mean(transformationFormulaAxis2$score / transformationFormulaAxis2$`-log10(q-value)`)

p <- VisualizeCNplot(scores_filt_call = scores_filt_call, GeneToVisualize = GeneToVisualize, chrArms = TRUE,
                      tabAnno_conf_Amp = tabAnno_conf_Amp, tabAnno_conf_Del = tabAnno_conf_Del,
                      FDR.thresh = FDR.thresh,
                      anno_Amp_GeneWidePeak = df_GR_ann_Amp_GeneWidePeak,
                      anno_Del_GeneWidePeak = df_GR_ann_Del_GeneWidePeak,
                      show.names = "all",
                      titleplot = paste0("TNBC", curSubtype, " subtype"))

p <- p + ggthemes::theme_base() + 
  labs(title = "", x = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

p


# Decompose the VisualizeCNplot function
scores_filt_call$Chromosome <- paste0("chr", scores_filt_call$Chromosome)

for (i in 1:9) {
  scores_filt_call[scores_filt_call$Chromosome == paste0("chr", i), "Chromosome"] <- paste0("chr0", i)
}
scores_filt_call <- cbind(scores_filt_call, ScoreSignif = rep(NA, nrow(scores_filt_call)))
scores_filt_call <- cbind(scores_filt_call, ChrPosLine = rep(0, nrow(scores_filt_call)))
scores_filt_call <- cbind(scores_filt_call, GeneAmp = rep("No", nrow(scores_filt_call)))
scores_filt_call <- cbind(scores_filt_call, GeneDel = rep("No", nrow(scores_filt_call)))

scores_filt_call$GeneAmp <- as.character(scores_filt_call$GeneAmp)
scores_filt_call$GeneDel <- as.character(scores_filt_call$GeneDel)

scores_filt_call_merged <- NULL

scores_filt_FDR_Amp <- scores_filt_call[scores_filt_call$`-log10(q-value)` > -log10(FDR.thresh) & scores_filt_call$Type %in% "Amp", ]
scores_filt_FDR_Del <- scores_filt_call[scores_filt_call$`-log10(q-value)` > -log10(FDR.thresh) & scores_filt_call$Type %in% "Del", ]

transformationFormulaAxis2 <- rbind(scores_filt_FDR_Amp, scores_filt_FDR_Del)
transformationFormulaAxis2_thresh <- mean(transformationFormulaAxis2$score / transformationFormulaAxis2$`-log10(q-value)`)

scoreThreshAmp <- min(scores_filt_FDR_Amp$score)
scoreThreshDel <- min(scores_filt_FDR_Del$score)

for (j in 1:23) {
  curChr <- paste0("chr0", j)
  if (j > 9) {
    curChr <- paste0("chr", j)
  }
  print(curChr)
  scores_filt_call_cur <- scores_filt_call[scores_filt_call$Chromosome %in% curChr, ]
  scores_filt_call_cur <- scores_filt_call_cur[order(scores_filt_call_cur$`Region Start [bp]`, decreasing = FALSE), ]
  
  if (nrow(scores_filt_call_cur) != 0) { 
    scores_filt_call_cur$`Region Start [bp]` == max(scores_filt_call_cur$`Region Start [bp]`)
    
  }
}

for (j in 1:23) {
  print(j)
}








scores_filt_call_cur <- scores_filt_call[scores_filt_call$Chromosome == 1,]
scores_filt_call_cur <- scores_filt_call[order(scores_filt_call$`Region Start [bp]`, decreasing = FALSE), ]
scores_filt_call_cur[scores_filt_call_cur$`Region Start [bp]` == max(scores_filt_call_cur$`Region Start [bp]`), "ChrPosLine"]

diamonds %>% group_by(clarity, cut) %>%
  ggplot(aes(x = clarity, y = price, group = cut, fill = cut)) + 
  geom_bar(stat = "identity")