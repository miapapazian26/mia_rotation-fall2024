library(cluster)  # load cluster package
library(ggplot2)  # load ggplot

# read in the data
data <- read.csv("GENENAMES_ecoli_expression_data.csv", header = TRUE)

# extract expression data (assumes columns 3 to the end contain the expression data)
expression_data <- data[, 3:ncol(data)]  
rownames(expression_data) <- data$Gene_name  # use gene names as row identifiers

# normalize the data
expression_data <- scale(expression_data)  

# calculate variance for each gene (row-wise)
gene_variance <- apply(expression_data, 1, var) #variance calculated for each gene row wise, higher variance = greater variability

# select the top 10 most variable genes
top_genes <- names(sort(gene_variance, decreasing = TRUE)[1:10]) #top ten genes w/ highest variance

# subset the data to only include the top 10 genes
expression_subset <- expression_data[top_genes, ] #filter

# perform CLARA clustering on the subset data
clara_result <- clara(expression_subset, k = 3)  # 3 clusters, adjust k as needed

# create a cluster data frame for the top genes
cluster_data <- data.frame(Gene_name = top_genes, Cluster = clara_result$clustering)

# merge the cluster labels back into the original data frame (only for top genes)
data <- merge(data, cluster_data, by = "Gene_name", all.x = TRUE)

# plot expression of the top 10 most variable genes
ggplot(data, aes(x = Gene_name, y = Rich_M9_ATC5_1.rep1, color = as.factor(Cluster))) +
  geom_point() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Gene Expression by Cluster (Top 10 Most Variable Genes)", color = "Cluster") +
  scale_x_discrete(limits = top_genes)  # limit the x-axis to the top genes


