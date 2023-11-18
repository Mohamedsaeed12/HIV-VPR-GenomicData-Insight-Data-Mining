#CLustering and Heatmap visualization 

#Data preparation
matching_rows <- normalized.expr$`Gene Symbol` %in% significant_gene_symbols
sig_gene_data <- normalized.expr[matching_rows, ]

#eucledian distance
numeric_data <- sig_gene_data[, -c(1, ncol(sig_gene_data))] # excluding the ID and Gene Symbol columns
dist_matrix <- dist(t(numeric_data))

#hierarchical clustering
gene_cluster <- hclust(dist_matrix)


library(gplots)
heatmap.2(as.matrix(numeric_data), 
          Rowv=as.dendrogram(gene_cluster), 
          scale="row", 
          trace="none", 
          margin=c(5,10), 
          cexRow=0.5, 
          cexCol=0.5, 
          col=colorRampPalette(c("blue", "white", "red"))(255))


dist_matrix <- dist(t(clustering_data))

gene_cluster <- hclust(dist_matrix)


summary(clustering_data)


clustering_data_numeric <- clustering_data[, -1] # Remove the first column (ID)

dist_matrix <- dist(t(clustering_data_numeric))
gene_cluster <- hclust(dist_matrix)
plot(gene_cluster, labels = FALSE, hang = -1, cex = 0.7)



# Assuming sig_gene_data is your data set and it's already preprocessed (normalized, NA values handled, etc.)

set.seed(123) # Setting seed for reproducibility

# Compute and plot wss (within-cluster sum of squares) for k = 1 to k = 10
wss <- sapply(1:10, function(k){
  kmeans(numeric_data, centers = k, nstart = 20)$tot.withinss
})

plot(1:10, wss, type = "b", xlab = "Number of clusters (k)", ylab = "Within groups sum of squares")



library(cluster)

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


#Kmeans clustering

set.seed(123) # for reproducibility
k <- 1 # for example, you can change this based on your data
kmeans_result <- kmeans(t(sig_gene_data[, -c(1, ncol(sig_gene_data))]), centers=k) # excluding the ID and Gene Symbol columns

ordered_data <- sig_gene_data[order(kmeans_result$cluster), ]


library(gplots)

heatmap.2(as.matrix(ordered_data[, -c(1, ncol(ordered_data))]), 
          scale="row", 
          trace="none", 
          margin=c(5,10), 
          cexRow=0.5, 
          cexCol=0.5, 
          col=colorRampPalette(c("blue", "white", "red"))(255),
          main=paste("K-means clustering with k =", k),
          labRow=ordered_data$`Gene Symbol`)

### DBSCAN 

install.packages("dbscan")
library(dbscan)

set.seed(123)  # For reproducibility
k <- 4  # Assuming minPts = 5
kNNdist <- kNNdistplot(numeric_data, k = k)
abline(h = 0.5, col = "red")  # Adjust the 0.5 based on the plot

dbscan_result <- dbscan(numeric_data, eps = 0.5, minPts = 4) 

# Assuming sig_gene_data has two columns for this example
plot(numeric_data, col=dbscan_result$cluster + 1L, pch=20, cex=1.5)
legend("topright", legend=unique(dbscan_result$cluster), col=1:length(unique(dbscan_result$cluster)), pch=20)


pca_result <- prcomp(numeric_data, scale. = TRUE)
plot(pca_result$x[,1:2], col=dbscan_result$cluster + 1L, pch=20, cex=1.5)


if (!requireNamespace("Rtsne", quietly = TRUE)) {
  install.packages("Rtsne")
}
library(Rtsne)
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

