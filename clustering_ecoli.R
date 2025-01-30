library(cluster)  # load cluster package
library(ggplot2)  # load ggplot

# read in the data
data <- read.csv("GENENAMES_ecoli_expression_data.csv", header = TRUE)

# extract expression data (assumes columns 3 to the end contain the expression data)
expression_data <- data[, 3:ncol(data)]  
rownames(expression_data) <- data$Gene_name  # use gene names as row identifiers

# normalize the data
expression_data <- scale(expression_data)  # or 

# apply log transformation to the data
expression_data_log <- log(expression_data + 1) #+1 to avoid log(0)

# calculate variance for each gene (row-wise)
gene_variance <- apply(expression_data, 1, var) # variance calculated for each gene row-wise, higher variance = greater variability

# select the top 10 most variable genes
top_genes <- names(sort(gene_variance, decreasing = TRUE)[1:10]) # top ten genes with highest variance

# subset the data to only include the top 10 genes
expression_subset <- expression_data[top_genes, ] # filter

# perform CLARA clustering on the subset data
clara_result <- clara(expression_subset, k = 3)  # 3 clusters, adjust k as needed
# or 
# clara_result <- clara(expression_data_log, k = 3)

# create a cluster data frame for the top genes
cluster_data <- data.frame(Gene_name = top_genes, Cluster = clara_result$clustering[match(top_genes, rownames(expression_subset))])

# merge the cluster labels for the top genes into the original data frame
data_top_genes <- merge(data[data$Gene_name %in% top_genes, ], cluster_data, by = "Gene_name", all.x = TRUE)

# plot expression of the top 10 most variable genes with their cluster labels
ggplot(data_top_genes, aes(x = Gene_name, y = log(Rich_M9_ATC5_1.rep1 + 1), color = as.factor(Cluster))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression by Cluster (Top 10 Most Variable Genes)", color = "Cluster") +
  scale_x_discrete(limits = top_genes)



