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


sig <- c('CSF','ApoE','TNF','APP','JAM','CCL')

# helper to extract grey-bar (total strength) for one object
get_grey_bar <- function(cc, sig, pattern="outgoing", slot.name="netP"){
  # make sure you've computed centrality
  if(length(cc@netP$centr)==0){
    cc <- netAnalysis_computeCentrality(cc, slot.name=slot.name)
  }
  centr <- cc@netP$centr
  # pull out the right measure and sum over cell types
  if(pattern=="outgoing"){
    tot <- sapply(sig, function(p) sum( centr[[p]]$outdeg ))
  } else if(pattern=="incoming"){
    tot <- sapply(sig, function(p) sum( centr[[p]]$indeg ))
  } else {
    tot <- sapply(sig, function(p) sum( centr[[p]]$outdeg + centr[[p]]$indeg ))
  }
  names(tot) <- sig
  tot
}

# apply to all of them
grey_list <- lapply(object.list, get_grey_bar, sig=sig)

# assemble into a single data.frame
grey_df <- do.call(cbind, grey_list)
rownames(grey_df) <- sig
colnames(grey_df) <- names(object.list)
View(grey_df)

# install.packages(c("tidyr","ggplot2","tibble"))
library(tibble)
library(tidyr)
library(ggplot2)

# reshape to long form
grey_long <- as.data.frame(grey_df) %>%
  rownames_to_column("pathway") %>%
  pivot_longer(
    cols      = -pathway,
    names_to  = "condition",
    values_to = "total"
  ) %>%
  # set factor levels to enforce plotting order:
  mutate(condition = factor(condition,
                            levels = c("WTB","A","GF","SPF","MBA","MBNA")))

# one plot per pathway
plots <- lapply(unique(grey_long$pathway), function(p) {
  dfp <- subset(grey_long, pathway == p)
  ggplot(dfp, aes(x = condition, y = total, fill = condition)) +
    geom_col() +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = 45, hjust = 1)
    ) +
    labs(
      title = p,
      x     = "Condition",
      y     = "Total outgoing strength"
    )
})

# display them
for (plt in plots) print(plt)



