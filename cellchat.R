yizhen <- LoadH5Seurat('G://yizhen_microglia_total_dataset_upd.h5Seurat')
head(yizhen)
yizhen <- RenameIdents(yizhen, 'Photoreceptor' = 'Rest','Photoreceptor_2' = 'Rest','18_Bipolar_Neuron' = 'Rest','Bipolar_Cone' = 'Rest','Microglia_Inflamm_C1'='Microglia','Microglia_IFN-Resp_C2' = 'Microglia', 'Microglia_resting_C0' = 'Microglia','Microglia_Proliferating_C3' = 'Microglia', '36_Melanocyte' = 'Melanocyte', '34_Melanocyte_2' = 'Melanocyte')
yizhen$ident_cellchat <- yizhen@active.ident

yizhen$ident_cellchat = droplevels(yizhen$ident_cellchat, exclude = setdiff(levels(yizhen$ident_cellchat),unique(yizhen$ident_cellchat)))
yizhen@active.ident = droplevels(yizhen@active.ident, exclude = setdiff(levels(yizhen@active.ident),unique(yizhen@active.ident)))

library(CellChat)
table(yizhen$Group)
yizhen_A <- subset(yizhen, subset = Group == 'A')
yizhen_A$ident_cellchat = droplevels(yizhen_A$ident_cellchat, exclude = setdiff(levels(yizhen_A$ident_cellchat),unique(yizhen_A$ident_cellchat)))
yizhen_A@active.ident = droplevels(yizhen_A@active.ident, exclude = setdiff(levels(yizhen_A@active.ident),unique(yizhen_A@active.ident)))
yizhen_KO <- subset(yizhen, subset = Group == 'KO')
yizhen_KO$ident_cellchat = droplevels(yizhen_KO$ident_cellchat, exclude = setdiff(levels(yizhen_KO$ident_cellchat),unique(yizhen_KO$ident_cellchat)))
yizhen_KO@active.ident = droplevels(yizhen_KO@active.ident, exclude = setdiff(levels(yizhen_KO@active.ident),unique(yizhen_KO@active.ident)))
yizhen_WTB1 <- subset(yizhen, subset = Group == 'B')
yizhen_WTB1$ident_cellchat = droplevels(yizhen_WTB1$ident_cellchat, exclude = setdiff(levels(yizhen_WTB1$ident_cellchat),unique(yizhen_WTB1$ident_cellchat)))
yizhen_WTB1@active.ident = droplevels(yizhen_WTB1@active.ident, exclude = setdiff(levels(yizhen_WTB1@active.ident),unique(yizhen_WTB1@active.ident)))
yizhen_WTB2 <- subset(yizhen, subset = Group == 'WT')
yizhen_WTB2$ident_cellchat = droplevels(yizhen_WTB2$ident_cellchat, exclude = setdiff(levels(yizhen_WTB2$ident_cellchat),unique(yizhen_WTB2$ident_cellchat)))
yizhen_WTB2@active.ident = droplevels(yizhen_WTB2@active.ident, exclude = setdiff(levels(yizhen_WTB2@active.ident),unique(yizhen_WTB2@active.ident)))
yizhen_WTB <- merge(yizhen_WTB1, y = yizhen_WTB2)
yizhen_GF <- subset(yizhen, subset = Group == 'GF')
yizhen_GF$ident_cellchat = droplevels(yizhen_GF$ident_cellchat, exclude = setdiff(levels(yizhen_GF$ident_cellchat),unique(yizhen_GF$ident_cellchat)))
yizhen_GF@active.ident = droplevels(yizhen_GF@active.ident, exclude = setdiff(levels(yizhen_GF@active.ident),unique(yizhen_GF@active.ident)))
yizhen_MBA <- subset(yizhen, subset = Group == 'MB-A')
yizhen_MBA$ident_cellchat = droplevels(yizhen_MBA$ident_cellchat, exclude = setdiff(levels(yizhen_MBA$ident_cellchat),unique(yizhen_MBA$ident_cellchat)))
yizhen_MBA@active.ident = droplevels(yizhen_MBA@active.ident, exclude = setdiff(levels(yizhen_MBA@active.ident),unique(yizhen_MBA@active.ident)))
yizhen_MBNA <- subset(yizhen, subset = Group == 'MB-NA')
yizhen_MBNA$ident_cellchat = droplevels(yizhen_MBNA$ident_cellchat, exclude = setdiff(levels(yizhen_MBNA$ident_cellchat),unique(yizhen_MBNA$ident_cellchat)))
yizhen_MBNA@active.ident = droplevels(yizhen_MBNA@active.ident, exclude = setdiff(levels(yizhen_MBNA@active.ident),unique(yizhen_MBNA@active.ident)))
yizhen_SPF <- subset(yizhen, subset = Group == 'SPF')
yizhen_SPF$ident_cellchat = droplevels(yizhen_SPF$ident_cellchat, exclude = setdiff(levels(yizhen_SPF$ident_cellchat),unique(yizhen_SPF$ident_cellchat)))
yizhen_SPF@active.ident = droplevels(yizhen_SPF@active.ident, exclude = setdiff(levels(yizhen_SPF@active.ident),unique(yizhen_SPF@active.ident)))

  cellchat <- createCellChat(object = yizhen_A , group.by = "ident_cellchat")
  CellChatDB <- CellChatDB.mouse
  CellChatDB.use <- CellChatDB
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- projectData(cellchat, PPI.mouse)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = 3)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

cellchat_yizhen_a <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_A.rds")

cellchat <- createCellChat(object = yizhen_KO , group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_KO <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_KO.rds")

cellchat <- createCellChat(object = yizhen_WTB , group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_WTB <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_WTB.rds")

cellchat <- createCellChat(object = yizhen_GF, group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_GF <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_GF.rds")

cellchat <- createCellChat(object = yizhen_MBA , group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_MBA <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_MBA.rds")

cellchat <- createCellChat(object = yizhen_MBNA , group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_MBNA <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_MBNA.rds")

cellchat <- createCellChat(object = yizhen_SPF , group.by = "ident_cellchat")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, raw.use = FALSE, population.size = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat_yizhen_SPF <- cellchat
saveRDS(cellchat, file = "G://YIZHEN/cellchat_SPF.rds")

ht1 <- netAnalysis_signalingRole_heatmap(cellchat_dogs_cui, pattern = "outgoing", height = 32)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat_dogs_cui, pattern = "incoming", height = 32)

ht <- netAnalysis_signalingRole_heatmap(cellchat_dogs_cui, pattern = "all", height = 32)
