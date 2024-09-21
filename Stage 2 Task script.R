url <- "https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv"
dataset <- read.csv(url)
library(gplots)

# Generate the heatmap
heatmap.2(as.matrix(dataset[, -1]),  # Assuming first column is gene names
          col = bluered(100),  # Diverging color palette
          trace = "none",
          dendrogram = "both",  # Cluster rows and columns
          margins = c(5, 6),
          scale = "row",
          main = "Heatmap of Gene Expression")
heatmap.2(as.matrix(dataset[,-1]), Rowv=T, Colv=F, dendrogram = "row", col = bluered(100), trace = "none", scale = "row") # CLUSTERED BY ROWS 

heatmap.2(as.matrix(dataset[,-1]), Rowv=F, Colv=T, dendrogram = "col", col = bluered(100), trace = "none", scale = "row") # CLUSTERED BY ROWS 

# Assuming 'dataset' is already loaded
# Set the first column name
colnames(dataset)[1] <- "Gene"

# Set the names for columns 2 to 11 as "1" to "10"
colnames(dataset)[2:11] <- as.character(1:10)

# Check the updated column names
print(colnames(dataset))


# Grouping the samples based on heat map cluster patterns
group1 <- data.frame(dataset[, c( 2, 3, 4, 5,6)])

group2 <- data.frame(dataset[, c( 7, 8, 9, 10,11)])

print(group2)

avg1 <- as.data.frame(rowMeans(group1))
avg2 <- as.data.frame(rowMeans(group2))


Foldchange <- log2(avg1) - log2(avg2)

p_values <- numeric(nrow(dataset))
for (i in 1:nrow(dataset))
{      
  g1 <- as.numeric(group1[i, ])  
  g2 <- as.numeric(group2[i, ])  
  
  t_test_result <- t.test(g1, g2, var.equal = FALSE)
  
  p_values[i] <- t_test_result$p.value
}

p_values <- as.data.frame(p_values)

print(p_values)

# Finding differentially regulated genes
genes <- data.frame(rownames(dataset))
genes <- data.frame(genes, p_values, Foldchange)
upregulated_genes <- subset(genes, Foldchange > 1.5 & p_values < 0.05)
downregulated_genes <- subset(genes, Foldchange < -1.5 & p_values < 0.05)

# Printing gene IDs 
downregulated_genes$rownames.dataset
upregulated_genes$rownames.dataset

install.packages("ggplot2")

library(ggplot2)


# Load necessary libraries
library(ggplot2)

# Create a bubble plot
ggplot(enrichment, aes(x = Fold.Enrichment, y = Pathway, size = nGenes)) +
  geom_point(alpha = 0.6) +  # Points with some transparency
  labs(title = "Enrichment Analysis Bubble Plot", 
       x = "Fold Enrichment", 
       y = "Pathway") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10))  # Customize y-axis text size






