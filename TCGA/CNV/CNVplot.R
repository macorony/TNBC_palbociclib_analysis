library(dplyr)
library(TCGAbiolinks)
library(GenomicRanges)
library(ggpubr)




load(file = "./data/TCGA/TCGA_TNBC_192samples_CNV.Rdata")
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
GeneToVisualize <- c(
  "PIK3CA","KDR","KIT","PDGFRA","FGFR1","BCL11A",
  "MLH1","APC","PTEN",
  "CASP8","PIK3CA","EGFR","JAK2",
  "ERRB3","KRAS","CTNNB1","MLH1","APC",
  "CSF1R","PTCH1","PTEN","RB1",
  "FGFR3","KDR","KIT","PDGFRA",
  "PTEN","FGFR2","BRCA2","RB1",
  "MAP2K1","BRCA1","PIK3CA",
  "FGFR2","MLH1","FGFR2","PTEN",
  "BRCA2",
  "CCNE1","AKT3","ATR","MLH1",
  "PIK3R1","ETV6","CREBBP","PALB2","NOTCH2",
  "NOTCH3","BRD4","MYC","B2M","HLA-A","HLA-B","HLA-C",
  "CCND2","FOXM1","MYB","PIM1","BAP1","MYCN","MCL1",
  "CHEK1","NFIB","CDKN2B","STK11","EZH2","SUZ12","SMARCD1",
  "PBRM1","DPF3","SMARCA4","ARID1B"
)


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
tabAnno_conf_Del <- readr::read_delim(paste0(folderCur, "/", fileCurDel), "\t", escape_double = FALSE, trim_ws = TRUE)
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

p <- VisualizeCNplot(scores_filt_call = scores_filt_call, GeneToVisualize = GeneToVisualize, chrArms = T,
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
    scores_filt_call_cur[scores_filt_call_cur$`Region Start [bp]` == max(scores_filt_call_cur$`Region Start [bp]`), "ChrPosLine"] <- curChr 
  }
  
  tabAnno_conf <- tabAnno_conf_Amp
  
  
  if (length(colnames(tabAnno_conf)) > 1) {
    curpos <- 1
    for (k in 2:(length(colnames(tabAnno_conf)) - 1)) {
      curSignif_Chr_q <- unlist(strsplit(as.character(colnames(tabAnno_conf)[k]), "q"))[1]
      
      curSignif_Chr_p <- unlist(strsplit(as.character(colnames(tabAnno_conf)[k]), "p"))[1]
      # Position p-arm
      if (curSignif_Chr_p == j) {
        if (length(intersect(scores_filt_call_cur$Type, "Amp")) == 1) {
          scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
          
          GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
          if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
            scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "GeneAmp"] <- 
              as.character(paste(intersect(GeneSignif_Amp_Del, GeneToVisualize), collapse = ";"))
          }
          curpos <- curpos + 1
          
        } 
      }
      # Position q-arm
      if (curSignif_Chr_q == j) {
        if (length(intersect(scores_filt_call_cur$Type, "Amp")) == 1) {
          scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
          
          GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
          
          if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
            scores_filt_call_cur[grep("Amp", scores_filt_call_cur$Type)[curpos], "GeneAmp"] <- 
              as.character(paste(intersect(GeneSignif_Amp_Del, GeneToVisualize), collapse = ";"))
          }
          curpos <- curpos + 1
          
        }
      }
    }
  }
  # Deletion file
  
  tabAnno_conf <- tabAnno_conf_Del
  
  if (length(colnames(tabAnno_conf)) > 1) {
    curpos <- 1
    for (k in 2:(length(colnames(tabAnno_conf)) - 1)) {
      curSignif_Chr_q <- unlist(strsplit(as.character(colnames(tabAnno_conf)[k]), "q"))[1]
      
      curSignif_Chr_p <- unlist(strsplit(as.character(colnames(tabAnno_conf)[k]), "p"))[1]
      # Position p-arm
      if (curSignif_Chr_p == j) {
        if (length(intersect(scores_filt_call_cur$Type, "Del")) == 1) {
          scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
          
          GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
          if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
            scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "GeneDel"] <- 
              as.character(paste(intersect(GeneSignif_Amp_Del, GeneToVisualize), collapse = ";"))
          }
          curpos <- curpos + 1
          
        } 
      }
      # Position q-arm
      if (curSignif_Chr_q == j) {
        if (length(intersect(scores_filt_call_cur$Type, "Del")) == 1) {
          scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "ScoreSignif"] <- colnames(tabAnno_conf)[k]
          
          GeneSignif_Amp_Del <- as.matrix(tabAnno_conf[, colnames(tabAnno_conf)[k]])
          
          if (length(intersect(GeneSignif_Amp_Del, GeneToVisualize)) != 0) {
            scores_filt_call_cur[grep("Del", scores_filt_call_cur$Type)[curpos], "GeneDel"] <- 
              as.character(paste(intersect(GeneSignif_Amp_Del, GeneToVisualize), collapse = ";"))
          }
          curpos <- curpos + 1
          
        }
      }
    }
  }
  scores_filt_call_merged <- rbind(scores_filt_call_merged, scores_filt_call_cur)
}

df <- subset(scores_filt_call_merged, select = c(
  "Type", "score", "Chromosome", "Region Start [bp]", "Region End [bp]", 
  "-log10(q-value)", "ScoreSignif", "ChrPosLine", "GeneAmp", "GeneDel"
))
df

colnames(df)[colnames(df) == "Region Start [bp]"] <- "start"
colnames(df)[colnames(df) == "Region End [bp]"] <- "end"

df[df$Type == "Del", "score"] <- -df[df$Type == "Del", "score"]

df <- cbind(Xpos = 1:nrow(df), df)

# Amp Genes

if (length(intersect(df$GeneAmp, GeneToVisualize)) == 0) {
  df_GR_ann_Amp <- 0
  length(df_GR_ann_Amp) <- 0
}

if (length(intersect(df$GeneAmp, GeneToVisualize)) != 0) {
  geneAmp_list <- NULL
  for (id in 1:nrow(df)) {
    curGene <- df$GeneAmp[id]
    
    if (curGene != "No") {
      geneAmp_list <- c(geneAmp_list, curGene)
    } 
  }
  geneAmp_list_split <- NULL
  
  for (id in 1:length(geneAmp_list)) {
    curGene <- geneAmp_list[id]
    
    if (curGene != "No") {
      geneAmp_list_split <- c(geneAmp_list_split, unlist(strsplit(curGene, ";")))
    }
  }
  geneAmp_list_split <- sort(geneAmp_list_split)
  
  df_GR_ann_Amp <- df_GR_ann_Amp_GeneWidePeak[df_GR_ann_Amp_GeneWidePeak$GeneSymbol %in% geneAmp_list_split, ]
  
  df_GR_ann_Amp <- df_GR_ann_Amp[order(df_GR_ann_Amp$score, decreasing = T), ]
  df_GR_ann_Amp <- df_GR_ann_Amp[!duplicated(df_GR_ann_Amp$GeneSymbol), ]
  
}


# Del Genes

if (length(intersect(df$GeneDel, GeneToVisualize)) == 0) {
  df_GR_ann_Del <- 0
  length(df_GR_ann_Del) <- 0
}

if (length(intersect(df$GeneDel, GeneToVisualize)) != 0) {
  geneDel_list <- NULL
  for (id in 1:nrow(df)) {
    curGene <- df$GeneDel[id]
    
    if (curGene != "No") {
      geneDel_list <- c(geneDel_list, curGene)
    } 
  }
  geneDel_list_split <- NULL
  
  for (id in 1:length(geneDel_list)) {
    curGene <- geneDel_list[id]
    
    if (curGene != "No") {
      geneDel_list_split <- c(geneDel_list_split, unlist(strsplit(curGene, ";")))
    }
  }
  geneDel_list_split <- unique(sort(geneDel_list_split))
  
  df_GR_ann_Del <- df_GR_ann_Del_GeneWidePeak[df_GR_ann_Del_GeneWidePeak$GeneSymbol %in% geneDel_list_split, ]
    
  df_GR_ann_Del <- df_GR_ann_Del[order(df_GR_ann_Del$score, decreasing = T), ]
  df_GR_ann_Del <- df_GR_ann_Del[!duplicated(df_GR_ann_Del$GeneSymbol), ]
  
}


df_mod <- df
df_mod$GeneAmp <- rep("No", nrow(df_mod))
df_mod$GeneDel <- rep("No", nrow(df_mod))

if (length(df_GR_ann_Del) != 0) {
  df_GR_ann_merged <- df_GR_ann_Del
}
if (length(df_GR_ann_Amp) != 0) {
  df_GR_ann_merged <- rbind(df_GR_ann_Amp, df_GR_ann_Del)
}

df_mod <- cbind(df_mod, 
                ChrStartEnd = paste0(gsub("chr", "", gsub("chr0", "", df_mod$Chromosome)), "_",
                                     df_mod$start, "_", df_mod$end)
                )

df_mod$ChrStartEnd <- as.character(df_mod$ChrStartEnd)

df_GR_ann_merged <- cbind(
  df_GR_ann_merged, 
  ChrStartEnd = paste0(
  df_GR_ann_merged$Chr, 
  "_",
  df_GR_ann_merged$start, 
  "_", 
  df_GR_ann_merged$end)
  )

df_GR_ann_merged$ChrStartEnd <- as.character(df_GR_ann_merged$ChrStartEnd)

show.names = "significant"
if (show.names == "significant") {
  df_GR_ann_merged <- df_GR_ann_merged[df_GR_ann_merged$`-log10(q-value)` > -log10(FDR.thresh), ]
  }

for (curSelIew in c("Amp", "Del")) {
  tabGRcur_ampDel <- df_GR_ann_merged[df_GR_ann_merged$Type %in% curSelIew, ]
  
  for (iw in 1:nrow(tabGRcur_ampDel)) {
    curgene <- tabGRcur_ampDel[iw, ]
    tabGrcur <- tabGRcur_ampDel[tabGRcur_ampDel$GeneSymbol %in% curgene$GeneSymbol, ]
    df_mod[df_mod$ChrStartEnd == tabGrcur$ChrStartEnd, paste0("Gene", curSelIew)] <- 
      curgene$GeneSymbol
  }
}
df <- df_mod
df$Chromosome <- factor(df$Chromosome, levels = unique(sort(df$Chromosome)))
require(ggrepel)

p1 <- ggplot() + geom_bar(aes(y = score, x = Xpos, col = Type), data = df, stat = "identity")
p1 <- p1 + scale_y_continuous(limits = c(-0.75, 1.5))
p1 <- p1 + geom_hline(yintercept = scoreThreshAmp, 
                      linetype = "dashed", 
                      color = "orange")
p1 <- p1 + geom_hline(yintercept = scoreThreshDel, 
                      linetype = "dashed",
                      color = "orange")

for (ichr in 1:23) {
  curPosChrsel <- df[df$Chromosome %in% paste0("chr0", ichr), ]
  if (ichr > 9) {
    curPosChrsel <- df[df$Chromosome %in% paste0("chr", ichr), ]
  }
  
  curPosChrsel <- curPosChrsel[curPosChrsel$ChrPosLine !=0, ]
  print(paste0("chr", ichr, " and ", curPosChrsel$Xpos))
  p1 <- p1 + geom_vline(
    xintercept = curPosChrsel$Xpos,
    linetype = "dashed",
    color = "black",
    alpha = 0.2
  )
}

chr.labels <- df %>% dplyr::filter(ChrPosLine != 0)
chr.labels
p1 <- p1 + scale_x_continuous(
  breaks = c(chr.labels$Xpos),
  limits = c(0, max(df$Xpos)),
  expand = c(0, 0), 
  labels = gsub("chr23", "chrX", c(paste0(chr.labels$Chromosome)))
)

p1 <- p1 + labs(
  title = paste0("Amplification and Deletion with Gistic2 in ", "TNBC all"),
  x = "Chromosomes",
  y = "G-scores"
)

df_Amp_gene <- df[df$Type %in% "Amp", ]
df_Amp_gene <- df_Amp_gene[df_Amp_gene$GeneAmp != "No", ]

df_Del_gene <- df[df$Type %in% "Del", ]
df_Del_gene <- df_Del_gene[df_Del_gene$GeneDel != "No", ]

if (nrow(df_Amp_gene) != 0) {
  p1 <- p1 + geom_point(
    data = df_Amp_gene, 
    shape = 1,
    color = "black",
    aes(x = Xpos, y = score)
  )
  
  p1 <- p1 + scale_shape(solid = FALSE)
  
  p1 <- p1 + geom_label_repel(
    data = df_Amp_gene, 
    aes(label = GeneAmp, x = Xpos, y = score), 
    nudge_y = 0.5,
    segment.size = 0.2,
    nudge_x = 0.5,
    box.padding = 0.35,
    point.padding = 0.5, 
    size = 3, 
    label.size = 0.1, 
    segment.color = "grey50"
  )
}

if (nrow(df_Del_gene) != 0) {
  p1 <- p1 + geom_point(
    data = df_Del_gene, 
    shape = 1,
    color = "black",
    aes(x = Xpos, y = score)
  )
  
  p1 <- p1 + geom_label_repel(
    data = df_Del_gene,
    aes(label = GeneDel, x = Xpos, y = score), 
    nudge_y = -0.75,
    nudge_x = -0.75,
    box.padding = 0.35,
    point.padding = 0.5,
    segment.size = 0.2,
    size = 3,
    label.size = 0.1,
    segment.color = "grey50"
  )
}



