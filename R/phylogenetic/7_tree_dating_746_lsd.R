
library(Rlsd2)

## date trees using lsd package


# wgs rate applied to whole tree- rate obtained from tree with WGS only
rate=0.000207203
all.result=lsd2(inputTree="output/phylogenetic/trees/all_rooted.746.nwk", inputDate ="output/phylogenetic/trees/all.dates.746.txt", outFile="output/phylogenetic/dated_trees/all_wgsrate_lsd_CI_746", seqLen=11925, ZscoreOutlier=5, givenRate=rate, confidenceInterval=1000)
#rate 0.000207203 [0.000207203; 0.000207203], tMRCA 1896.57 [1885.48; 1906.56]

