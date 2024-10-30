#Load adult human Seurat object
library(reshape2)
library(tidyr)
library(tidyverse)
library(dplyr)

avg_expr <- AverageExpression(adult, features = c('IDO1','IDO2','KYNU','KMO','QPRT','KYAT','LAT1'))
avg_expr <- AverageExpression(adult, features = c('SLC7A5'))
avg_expr_matrix <- avg_expr$RNA
avg_expr_long <- melt(as.matrix(avg_expr_matrix))
colnames(avg_expr_long) <- c("Gene", "Group", "Expression")

ggplot(avg_expr_long, aes(x = Group, y = Gene, fill = Expression)) +
geom_tile() +
scale_fill_gradientn(colours = c("blue", "white", "red")) + theme_minimal() +
# Custom color gradient
theme(axis.text.x = element_text(angle = 90, hjust = 1)) + # Adjust x-axis text
labs(x = "Group", y = "Gene", fill = "Expression") + coord_fixed()
