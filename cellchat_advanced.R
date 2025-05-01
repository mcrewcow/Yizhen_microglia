pathway = 'TNF'
netAnalysis_contribution(cellchat_yizhen_WTB, signaling = pathway)
netAnalysis_contribution(cellchat_yizhen_a, signaling = pathway)
netAnalysis_contribution(cellchat_yizhen_GF, signaling = pathway)
netAnalysis_contribution(cellchat_yizhen_SPF, signaling = pathway)
netAnalysis_contribution(cellchat_yizhen_MBA, signaling = pathway)
netAnalysis_contribution(cellchat_yizhen_MBNA, signaling = pathway)

contrib_list <- lapply(names(object.list), function(nm) {
  res <- netAnalysis_contribution(
    object      = object.list[[nm]],
    signaling   = pathway,
    return.data = TRUE
  )
  df <- res$LR.contribution
  df$Condition <- nm
  df
})

# combine into one big data.frame
library(dplyr)
contrib_df <- bind_rows(contrib_list)

# inspect
print(contrib_df)

contrib_df$Condition <- factor(
  contrib_df$Condition,
  levels = c("WTB","A","GF","SPF","MBA","MBNA")
)

# 2) (optional) enforce a specific order of L-R pairs along the x-axis;
#    by default this will be alphabetical, but you can customize:
pair_order <- unique(contrib_df$name)  
contrib_df$name <- factor(contrib_df$name, levels = pair_order)

# 3) plot
ggplot(contrib_df, aes(x = name, y = contribution, fill = Condition)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_brewer(palette = "Set2") +     # you can pick any palette you like
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  ) +
  labs(
    x     = "Ligand–Receptor Pair",
    y     = "Relative Contribution",
    fill  = "Condition",
    title = paste0("TNF pathway contributions")
  )





