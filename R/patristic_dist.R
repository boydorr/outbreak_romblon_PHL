## patristic distances for sequenced romblon cases
library(ape)
library(adephylo)
library(reshape2)
library(dendextend)
library(RColorBrewer)

# import tree
romblon.tree=read.tree("output/phylogenetic_distances/romblonSeq_24.nwk")

# patristic distances
pat=distTips(romblon.tree, tips = "all", method ="patristic",useC = TRUE)
pat.mat=as.matrix(pat)

write.table(pat.mat, file="output/phylogenetic_distances/romblonSeq_24_patristicDist_matrix.txt")

# Plot frequency distribution of Hamming distances
hist(pat.mat, 
     breaks = seq(0, max(pat.mat) + 0.0005, by=0.0005),
     col = "cyan",
     xlab = "Hamming distance", main = "Distribution of frequencies")
dev.off()

heatmap(pat.mat, col = brewer.pal(9, "Blues")) 
pat.melt <- melt(pat.mat)

# Perform hierarchical clustering
hc <- hclust(as.dist(pat.mat))

# Get the order of tips from the dendrogram
dend_order <- order.dendrogram(as.dendrogram(hc))

# Reorder the matrix based on the dendrogram order
pat.mat.ordered <- pat.mat[dend_order, dend_order]

# Melt the ordered matrix for ggplot
pat.melt <- melt(pat.mat.ordered)

# Create the heatmap with ggplot
ggplot(pat.melt, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) + 
  theme_minimal() +
  labs(x = "Tip 1", y = "Tip 2", fill = "Distance",
       title = "Patristic Distance Heatmap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_colorbar(title = "Distance"))


###same for whole tree (n=518)
whole.tree=read.tree("output/phylogenetic/trees/all_rooted.nwk")

pat2=distTips(whole.tree, tips = "all", method ="patristic",useC = TRUE)
pat.mat2=as.matrix(pat2)

write.table(pat.mat2, file="output/phylogenetic_distances/contextualTree_patristicDist_matrix.txt")

# Plot frequency distribution of Hamming distances
hist(pat.mat2, 
     breaks = seq(0, max(pat.mat2) + 0.0005, by=0.0005),
     col = "cyan",
     xlab = "Hamming distance", main = "Distribution of frequencies")
dev.off()

heatmap(pat.mat2, col = brewer.pal(9, "Blues")) 
pat.melt2 <- melt(pat.mat2)

# Perform hierarchical clustering
hc2 <- hclust(as.dist(pat.mat2))

# Get the order of tips from the dendrogram
dend_order2 <- order.dendrogram(as.dendrogram(hc2))

# Reorder the matrix based on the dendrogram order
pat.mat.ordered2 <- pat.mat[dend_order2, dend_order2]

# Melt the ordered matrix for ggplot
pat.melt2 <- melt(pat.mat.ordered2)

# Create the heatmap with ggplot
ggplot(pat.melt2, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "Blues")) + 
  theme_minimal() +
  labs(x = "Tip 1", y = "Tip 2", fill = "Distance",
       title = "Patristic Distance Heatmap") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=1),axis.text.y = element_text(size=1)) +
  guides(fill = guide_colorbar(title = "Distance"))



