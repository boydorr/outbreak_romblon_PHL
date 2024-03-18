library(adegenet)
library(here)

## Cluster analysis of genetic data, and lineage assignments

pat_matrix <- as.matrix(read.table(here("output/phylogenetic_distances/romblonSeq_24_patristicDist_matrix.txt"), check.names = FALSE))

# Plot frequency distribution of Hamming distances
jpeg(here("output/figures/Hamming_distance_frequencies.jpeg"))
hist(pat_matrix, 
     breaks = seq(0, max(pat_matrix) + 0.0005, by=0.0005),
     col = "cyan",
     xlab = "Hamming distance", main = "Distribution of frequencies")
dev.off()

# GENETIC CLUSTERS- ADEGENET
# The resulting cut-off of 0.0004 was selected by comparing resulting clusters with
# phylogenetic trees .
clust_all <- gengraph(pat_matrix, show.graph=F, col.pal=rainbow, truenames=T, cex=0.3, cutoff=0.0004)
write.csv(clust_all$clust$membership, here("output/phylogenetic_distances/0.0004clust.csv"))
