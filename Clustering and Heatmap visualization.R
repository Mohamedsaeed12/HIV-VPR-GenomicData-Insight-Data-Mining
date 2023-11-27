# Clustering and Heatmap visualization 

# Load necessary libraries
library(gplots)
library(cluster)
library(dbscan)
library(Rtsne)

# Data preparation
matching_rows <- normalized.expr$`Gene Symbol` %in% significant_gene_symbols
sig_gene_data <- normalized.expr[matching_rows, ]

# Euclidean distance
numeric_data <- sig_gene_data[, -c(1, ncol(sig_gene_data))] # excluding the ID and Gene Symbol columns
dist_matrix <- dist(t(numeric_data))

# Hierarchical clustering
gene_cluster <- hclust(dist_matrix)

# Assuming that numeric_data has 4 columns, corresponding to the conditions and donors as mentioned.
# Rename the columns as specified.
colnames(numeric_data) <- c("Vpr donor 1","Zs-Mac donor 1", "Vpr donor 2", "Zs-Mac donor 2")

# Heatmap visualization
heatmap.2(as.matrix(numeric_data), 
          Rowv=as.dendrogram(gene_cluster), 
          scale="row", 
          trace="none", 
          margin=c(10,3),  # Adjusted margin to accommodate row and column names
          cexRow=0.5, 
          cexCol=1,  # Increased to make column names more readable
          col=colorRampPalette(c("blue", "white", "red"))(255),
          main="HC Heatmap",
          labCol=colnames(numeric_data))  # Ensure column names are added as labels


# Assuming sig_gene_data is your data set and it's already preprocessed (normalized, NA values handled, etc.)
set.seed(123) # Setting seed for reproducibility

# Compute and plot wss (within-cluster sum of squares) for k = 1 to k = 10
wss <- sapply(1:10, function(k) {
  kmeans(numeric_data, centers = k, nstart = 20)$tot.withinss
})

plot(1:10, wss, type = "b", xlab = "Number of clusters (k)", ylab = "Within groups sum of squares")

# Compute average silhouette for different numbers of clusters
sil_width <- sapply(2:10, function(k) {
  model <- kmeans(numeric_data, centers = k, nstart = 20)
  silhouette_score <- silhouette(model$cluster, dist(numeric_data))
  mean(silhouette_score[, 3])
})

# Find the optimal number of clusters with the maximum average silhouette width
optimal_k <- which.max(sil_width)

# Plot the silhouette width for each k
plot(2:10, sil_width, type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters k", ylab = "Average silhouette width")
abline(v = optimal_k, lty = 2)

# Print the optimal k
cat("The optimal number of clusters k is:", optimal_k, "\n")

# Kmeans clustering
set.seed(123) # for reproducibility
k <- optimal_k # set to the optimal number of clusters
kmeans_result <- kmeans(t(numeric_data), centers=k) # transposed to cluster samples instead of genes

# Ordering data based on kmeans clusters
ordered_data <- sig_gene_data[order(kmeans_result$cluster), ]

colnames(ordered_data) <- c("ID","Vpr donor 1","Zs-Mac donor 1", "Vpr donor 2", "Zs-Mac donor 2")
# Heatmap visualization for K-means clustering
heatmap.2(as.matrix(ordered_data[, -c(1, ncol(ordered_data))]), 
          scale="row", 
          trace="none", 
          margin=c(5,3), 
          cexRow=0.5, 
          cexCol=0.5, 
          col=colorRampPalette(c("blue", "white", "red"))(255),
          main=paste("KC with k =", k),
          labRow=ordered_data$`Gene Symbol`) # Check if this column exists for labels

# DBSCAN clustering
set.seed(123)  # For reproducibility
k <- 4  # Assuming minPts = 5 for kNNdistplot
kNNdist <- kNNdistplot(numeric_data, k = k)
abline(h = 0.5, col = "red")  # Adjust the 0.5 based on the plot

# You might need to choose a different eps based on the kNNdist plot
dbscan_result <- dbscan(numeric_data, eps = 0.5, minPts = k)

# PCA for visualization
pca_result <- prcomp(numeric_data, scale. = TRUE)
plot(pca_result$x[,1:2], col=dbscan_result$cluster + 1L, pch=20, cex=1.5)

# t-SNE for visualization
set.seed(123)  # For reproducibility
tsne_result <- Rtsne(numeric_data, dims = 2, perplexity = 30, theta = 0.5, max_iter = 1000)
plot(tsne_result$Y, col=dbscan_result$cluster + 1L, asp=1, pch=20, cex=1.5)
legend("topright", legend=unique(dbscan_result$cluster), col=1:length(unique(dbscan_result$cluster)), pch=20)

# Define a color for each cluster and a distinct color for noise points
cluster_colors <- rainbow(length(unique(dbscan_result$cluster)))
noise_color <- "black"  # Color for noise points
colors <- ifelse(dbscan_result$cluster == -1, noise_color, cluster_colors[dbscan_result$cluster])

plot(tsne_result$Y, col=colors, asp=1, pch=20, cex=1.5)
legend("topright", legend=c("Noise", unique(dbscan_result$cluster)), col=c(noise_color, cluster_colors), pch=20)

