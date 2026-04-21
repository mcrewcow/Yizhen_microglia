yizhen_mg <- LoadH5Seurat('C://Bioinf/Yizhen/escape_GSEA_UPD_v3_MICROGLIA_ANNO.h5Seurat')


ES2 <- data.frame(yizhen_mg[[]], Idents(yizhen_mg))
colnames(ES2)[ncol(ES2)] <- "cluster"

head(ES2)
ES2$Group <- as.character(ES2$Group)
ES2$Group[ES2$Group %in% c("WT", "B")] <- "WT"
ES2$Group <- factor(ES2$Group)

ES2$Group2 <- as.character(ES2$Group)

ES2$Group2 <- dplyr::case_when(
  ES2$Group2 %in% c("WT", "B") ~ "WT",
  ES2$Group2 %in% c("MB-NA", "MBNA") ~ "MBNA",
  ES2$Group2 %in% c("MB-A", "MBA") ~ "MBA",
  TRUE ~ ES2$Group2
)

ES2$Group2 <- factor(
  ES2$Group2,
  levels = rev(c("WT", "A", "MBNA", "MBA", "SPF", "GF"))
)

p <- ridgeEnrichment(
  ES2,
  gene.set = "GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY",
  group = "Group2",
  add.rug = FALSE
)
p

cut_first_third_left <- 750
cut_first_third_right <- 2400

s <- ES2$GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY

ES2$peak_manual <- NA_character_
ES2$peak_manual[s <= cut_first_third_left] <- "first_peak"
ES2$peak_manual[s >= cut_first_third_right] <- "third_peak"

yizhen_mg$peak_manual <- ES2$peak_manual
table(yizhen_mg$peak_manual, useNA = "ifany")

yizhen_mg$peak_manual <- factor(yizhen_mg$peak_manual, levels = c("first_peak", "third_peak"))
Idents(yizhen_mg) <- yizhen_mg$peak_manual
deg_first_vs_third <- FindMarkers(
  yizhen_mg,
  ident.1 = "third_peak",
  ident.2 = "first_peak",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(deg_first_vs_third)
table(yizhen_mg$peak_manual, useNA = "ifany")
mg_peak <- subset(
  yizhen_mg,
  subset = peak_manual %in% c("first_peak", "third_peak")
)

table(mg_peak$peak_manual, useNA = "ifany")
mg_peak$peak_manual <- factor(
  mg_peak$peak_manual,
  levels = c("third_peak", "first_peak")
)
DotPlot(mg_peak, features = c( 'Trem2', 
                               'C1qa', 'C1qb', 'Ctss', 'Ctsd', 'Tyrobp', 'Apoe', 'Csf1r','Treml4', 'Ifitm6', 'Trem3', 'Clec4e', 'Cd274', 'Itga4', 'Trem1'), group.by = 'peak_manual', assay = 'RNA')

p <- ridgeEnrichment(
  ES2,
  gene.set = "GOBP_TRYPTOPHAN_TRANSPORT",
  group = "Group2",
  add.rug = FALSE
)
p

cut_first_third_left <- -3300
cut_first_third_right <- 400

s <- ES2$GOBP_TRYPTOPHAN_TRANSPORT

ES2$peak_manual <- NA_character_
ES2$peak_manual[s <= cut_first_third_left] <- "first_peak"
ES2$peak_manual[s >= cut_first_third_right] <- "third_peak"

yizhen_mg$peak_manual <- ES2$peak_manual
table(yizhen_mg$peak_manual, useNA = "ifany")

yizhen_mg$peak_manual <- factor(yizhen_mg$peak_manual, levels = c("first_peak", "third_peak"))
Idents(yizhen_mg) <- yizhen_mg$peak_manual
deg_first_vs_third <- FindMarkers(
  yizhen_mg,
  ident.1 = "third_peak",
  ident.2 = "first_peak",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

head(deg_first_vs_third)
table(yizhen_mg$peak_manual, useNA = "ifany")
mg_peak <- subset(
  yizhen_mg,
  subset = peak_manual %in% c("first_peak", "third_peak")
)

table(mg_peak$peak_manual, useNA = "ifany")
mg_peak$peak_manual <- factor(
  mg_peak$peak_manual,
  levels = c("third_peak", "first_peak")
)
DotPlot(mg_peak, features = c( 'Trem2', 
                               'C1qa', 'C1qb', 'Ctss','C1qc', 'Aif1',  'Tyrobp', 'Apoe', 'Csf1r',
                               'Il1a', 'Jun', 'Fosb', 'Ccl4', 'Ccl3', 'Tnf', 'Junb', 'Nktr', 'Socs3', 'Jund',
                               'Ccl2',  'Ctsb', 'Cx3cr1','Ctsd'), group.by = 'peak_manual', assay = 'RNA')
