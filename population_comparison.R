avg_expr <- AverageExpression(yizhen, features = c('Ido1','Ido2','Kynu','Kmo','Qprt','Kyat','Lat1'), group.by = c('predicted.id','Group'))
library(ggplot2)
library(reshape2)
avg_expr <- AverageExpression(yizhen_mg, features = c('Ido1','Ido2','Kynu','Kmo','Qprt','Kyat','Lat1'), group.by = c('microglia_state','Group'))

avg_expr_matrix <- avg_expr$RNA

# Convert matrix to long format for ggplot2
avg_expr_long <- melt(as.matrix(avg_expr_matrix))

# Rename the columns for better understanding
colnames(avg_expr_long) <- c("Gene", "Group", "Expression")

# Create the heatmap using ggplot2
ggplot(avg_expr_long, aes(x = Group, y = Gene, fill = Expression)) +
  geom_tile() +
  scale_fill_gradientn(colours = c("blue", "white", "red")) + # Custom color gradient
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Adjust x-axis text
  labs(x = "Group", y = "Gene", fill = "Expression") + coord_fixed()

# Calculate the 25th and 75th percentiles
percentile_25 <- quantile(yizhen_mg$GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY, 0.25) #GOBP_TRYPTOPHAN_TRANSPORT
percentile_75 <- quantile(yizhen_mg$GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY, 0.75)

# Create a new column with 'positive' or 'negative' labels
yizhen_mg$Tryptophan_Transport_Status <- ifelse(
  yizhen_mg$GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY >= percentile_75, "positive", 
  ifelse(yizhen_mg$GOBP_IMMUNE_RESPONSE_REGULATING_SIGNALING_PATHWAY <= percentile_25, "negative", NA)
)

# Check the results
table(yizhen_mg$Tryptophan_Transport_Status)
yizhen_mgpos <-  subset(yizhen_mg, subset = Tryptophan_Transport_Status == 'positive')
yizhen_mgneg <-  subset(yizhen_mg, subset = Tryptophan_Transport_Status == 'negative')
yizhen_mgposneg <- merge(x = yizhen_mgpos, y = yizhen_mgneg)
table(yizhen_mgposneg$Tryptophan_Transport_Status)
yizhen_mgposneg <- SetIdent(yizhen_mgposneg, value = 'Tryptophan_Transport_Status')
markers <- FindAllMarkers(yizhen_mgposneg, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
markers  %>%
group_by(cluster) %>%
slice_max(n = 10, order_by = avg_log2FC)
markers  %>%
group_by(cluster) %>%
top_n(n = 20, wt = avg_log2FC) -> top10
DoHeatmap(yizhen_mgposneg, features = top10$gene) + NoLegend()

DotPlot(yizhen_mgposneg, features = c('Slc3a2','Slc7a8','Egr1','Ccl4','Atf3','Fosb','Ccl3','Il1a','Jun','Fos','Ccl5'))

DotPlot(yizhen_mgposneg, features = c('Atf3','Klf6','Ccl4','Btg2','Klf2','Egr1','Jun','Nfkbiz','Nfkbia','Ccl3','Rgs10','Gng10','Fau'))
