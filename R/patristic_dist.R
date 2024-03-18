## patristic distances for sequenced romblon cases
library(ape)
library(adephylo)

# import tree
romblon.tree=read.tree("processed_data/romblon_sequences/romblonSeq_24.nwk")

# patristic distances
pat=distTips(romblon.tree, tips = "all", method ="patristic",useC = TRUE)
pat.mat=as.matrix(pat)

write.table(pat.mat, file="processed_data/romblon_sequences/romblonSeq_24_patristicDist_matrix.txt")

# Plot frequency distribution of Hamming distances
hist(pat_matrix, 
     breaks = seq(0, max(pat_matrix) + 0.0005, by=0.0005),
     col = "cyan",
     xlab = "Hamming distance", main = "Distribution of frequencies")
dev.off()
