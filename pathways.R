library(Seurat)
library(SeuratDisk)

head(yizhen)
yizhen <- SetIdent(yizhen, value = 'orig.ident')
table(yizhen$orig.ident)
yizhen <- RenameIdents(yizhen, 'WT' = 'WTB', 'B' = 'WTB')
yizhen$background <- yizhen@active.ident

library(escape)

gene.sets1 <- getGeneSets(library = "C5", gene.sets = c("GOBP_MICROGLIA_DIFFERENTIATION",
'GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE','GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY',
'GOBP_MICROGLIAL_CELL_MIGRATION','GOBP_MICROGLIAL_CELL_PROLIFERATION','GOBP_REGULATION_OF_MICROGLIAL_CELL_ACTIVATION',
'GOBP_MACROPHAGE_ACTIVATION','GOBP_MACROPHAGE_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
'GOBP_MACROPHAGE_APOPTOTIC_PROCESS',
'GOBP_MACROPHAGE_CHEMOTAXIS',
'GOBP_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION',
'GOBP_MACROPHAGE_COLONY_STIMULATING_FACTOR_SIGNALING_PATHWAY',
'GOBP_MACROPHAGE_CYTOKINE_PRODUCTION',
'GOBP_MACROPHAGE_DIFFERENTIATION',
'GOBP_MACROPHAGE_FUSION',
'GOBP_MACROPHAGE_INFLAMMATORY_PROTEIN_1_ALPHA_PRODUCTION',
'GOBP_MACROPHAGE_MIGRATION',
'GOBP_MACROPHAGE_PROLIFERATION',
'GOBP_REGULATION_OF_MACROPHAGE_ACTIVATION',
'GOBP_REGULATION_OF_MACROPHAGE_APOPTOTIC_PROCESS',
'GOBP_REGULATION_OF_MACROPHAGE_CHEMOTAXIS',
'GOBP_REGULATION_OF_MACROPHAGE_DERIVED_FOAM_CELL_DIFFERENTIATION',
'GOBP_REGULATION_OF_MACROPHAGE_DIFFERENTIATION',
'GOBP_REGULATION_OF_MACROPHAGE_FUSION',
'GOBP_REGULATION_OF_MACROPHAGE_MIGRATION',
'GOBP_REGULATION_OF_MACROPHAGE_PROLIFERATION',
'HP_ABNORMAL_MACROPHAGE_MORPHOLOGY',
'GOBP_RESPONSE_TO_MACROPHAGE_COLONY_STIMULATING_FACTOR',
'GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_MACROPHAGE_COLONY_STIMULATING_FACTOR_STIMULUS',
'GOBP_RESPONSE_TO_MACROPHAGE_COLONY_STIMULATING_FACTOR',
'GOBP_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION',
'GOBP_POSITIVE_REGULATION_OF_GRANULOCYTE_MACROPHAGE_COLONY_STIMULATING_FACTOR_PRODUCTION',
'GOBP_PHAGOCYTOSIS',
'GOBP_PHAGOCYTOSIS_ENGULFMENT',
'GOBP_PHAGOCYTOSIS_RECOGNITION',
'GOBP_T_CELL_ACTIVATION',
'GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
'GOBP_T_CELL_ACTIVATION_VIA_T_CELL_RECEPTOR_CONTACT_WITH_ANTIGEN_BOUND_TO_MHC_MOLECULE_ON_ANTIGEN_PRESENTING_CELL',
'GOBP_T_CELL_ANTIGEN_PROCESSING_AND_PRESENTATION',
'GOBP_T_CELL_APOPTOTIC_PROCESS',
'GOBP_T_CELL_CHEMOTAXIS',
'GOBP_T_CELL_CYTOKINE_PRODUCTION',
'GOBP_T_CELL_DIFFERENTIATION'),species = 'Mus musculus')

ES <- enrichIt(obj = yizhen,
gene.sets = gene.sets1,
groups = 1000, cores = 1)

yizhen <- AddMetaData(yizhen, ES)
ES <- data.frame(yizhen[[]], Idents(yizhen))
colnames(ES)[ncol(ES)] <- "cluster"
ridgeEnrichment(ES, gene.set = "MP", group = 'background', add.rug = TRUE)


ES_macrophage <- subset(ES, subset = ident_v1 == '31_Macrophage')
ridgeEnrichment(ES_macrophage, gene.set = "GOBP_MACROPHAGE_ACTIVATION", group = 'background', add.rug = TRUE) + facet_wrap(~ident_v1)

yizhen <- RenameIdents(yizhen, 'Photoreceptor' = 'Rest','Photoreceptor_2' = 'Rest','18_Bipolar_Neuron' = 'Rest','Bipolar_Cone' = 'Rest','Microglia_Inflamm_C1'='Microglia','Microglia_IFN-Resp_C2' = 'Microglia', 'Microglia_resting_C0' = 'Microglia','Microglia_Proliferating_C3' = 'Microglia', '36_Melanocyte' = 'Melanocyte', '34_Melanocyte_2' = 'Melanocyte')
yizhen$ident_v2 <- yizhen@active.ident
yizhen_mg <- subset(yizhen, idents = c('Microglia','31_Macrophage'))

ES_mg <- subset(ES, subset = ident_v2 == c('Microglia','31_Macrophage'))
ridgeEnrichment(ES_mg, gene.set = "GOBP_MICROGLIAL_CELL_MEDIATED_CYTOTOXICITY", group = 'background', add.rug = TRUE) + facet_wrap(~ident_v2)

ES_T <- subset(ES, subset = ident_v2 == 'T_Cell')

ridgeEnrichment(ES_T, gene.set = "GOBP_T_CELL_ACTIVATION", group = 'background', add.rug = TRUE) + facet_wrap(~ident_v2)

