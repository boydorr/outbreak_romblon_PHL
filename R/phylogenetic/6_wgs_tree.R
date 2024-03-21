library(seqinr)
library(ggtree)
#devtools::install_github("tothuhien/Rlsd2")
library(Biostrings)
library(stringr)
library(BactDating)
library(ape)
library(castor)

# prep file that only contains WGS (i.e. >90% genome coverage)
# count lengths of all seq in files (without gaps or Ns)
string=readDNAStringSet("processed_data/concatenated_alignment/ph_all_merged_by_id_581.fasta")
lengths=str_count(as.character(string), "A|T|C|G")
wgs.names=names(string)[which(lengths>max(lengths)*0.9)]
wgs.names=data.frame(Accession=wgs.names)
write.table(wgs.names, file = "processed_data/wgs_alignment/wgs.names.txt", sep = "", row.names=F, quote=F)

## wgs associated metadata
all_phl.annot=read.csv("processed_data/concatenated_alignment/ph_metadata_merged_by_id_581.csv")
wgs.metadata=all_phl.annot[which(all_phl.annot$Sample.ID %in% wgs.names$Accession),]
write.csv(wgs.metadata, file = "processed_data/wgs_alignment/wgs.metadata.csv", row.names=F)

## root the tree with all philippines cases
all_tree=read.tree("processed_data/trees/ph_all_581_ft.nwk")
# create dates file for lsd2
# important - lsd2 in R requires first line to be number of sequences
all.dates <- data.frame(
  Sample.ID = all_phl.annot$Sample.ID,
  decimal_date = all_phl.annot$decimal_date
)
seqs=nrow(all_phl.annot)
# Create a new row dataframe with the value of seqs
new_row <- data.frame(Sample.ID = seqs, decimal_date = "")
# Use rbind to bind the new row to the top of wgs.dates
all.dates <- rbind(new_row,all.dates)
write.table(all.dates, file = "processed_data/trees/all.dates.txt", row.names=F, sep="\t", quote = F,col.names = F)

# re-root the tree based on rtt regression
ordered.dates=all_phl.annot$decimal_date[match(all_tree$tip.label,all_phl.annot$Sample.ID)]
rooted=initRoot(all_tree,as.numeric(ordered.dates))
ordered.dates2=all_phl.annot$decimal_date[match(rooted$tip.label,all_phl.annot$Sample.ID)]
rooted$edge.length=rooted$edge.length*11925
r=roottotip(rooted,as.numeric(ordered.dates2))
rooted$edge.length=rooted$edge.length/11925
write.tree(rooted,"processed_data/trees/all_rooted.nwk")
# outside of R- use gotree prune to extract WGS only tree 

# create dates file for lsd2
# important - lsd2 in R requires first line to be number of sequences
wgs.dates <- data.frame(
  Sample.ID = wgs.metadata$Sample.ID,
  decimal_date = wgs.metadata$decimal_date
)
seqs=nrow(wgs.dates)
# Create a new row dataframe with the value of seqs
new_row <- data.frame(Sample.ID = seqs, decimal_date = "")
# Use rbind to bind the new row to the top of wgs.dates
wgs.dates <- rbind(new_row, wgs.dates)
write.table(wgs.dates, file = "processed_data/wgs_alignment/wgs.dates.txt", row.names=F, sep="\t", quote = F,col.names = F)

# create a clockor2 file
wgs.dates=cbind(tip=wgs.metadata$Sample.ID,date=wgs.metadata$decimal_date, group=wgs.metadata$Province)
write.csv(wgs.dates, file = "processed_data/wgs_alignment/wgs.clockor2.csv", row.names=F, quote = F)

## import the wgs tree and root according to best temporal root 
### trees
wgs.tree=read.tree("processed_data/wgs_alignment/ph_wgs_ft.nwk")
ordered.dates=wgs.metadata$decimal_date[match(wgs.tree$tip.label,wgs.metadata$Sample.ID)]
rooted=initRoot(wgs.tree,as.numeric(ordered.dates))
ordered.dates2=wgs.metadata$decimal_date[match(rooted$tip.label,wgs.metadata$Sample.ID)]
rooted$edge.length=rooted$edge.length*11925
r=roottotip(rooted,as.numeric(ordered.dates2))
rooted$edge.length=rooted$edge.length/11925
write.tree(rooted,"processed_data/wgs_alignment/wgs_rooted.nwk")


result=lsd2(inputTree="processed_data/wgs_alignment/wgs_rooted.nwk", inputDate ="processed_data/wgs_alignment/wgs.dates.txt", outFile="processed_data/wgs_alignment//wgstest", seqLen=11925, ZscoreOutlier=3)

result2=lsd2(inputTree="processed_data/wgs_alignment/wgs_only_prerooted_tree.nwk", inputDate ="processed_data/wgs_alignment/wgs.dates.txt", outFile="processed_data/wgs_alignment/wgs_prerooted", seqLen=11925, ZscoreOutlier=3)

# wgs rate applied to whole tree
rate=0.000203983
# lsd2 -i data/TempestRooted_RABV_canine.nwk -d data/canine_lsd2_dates.txt -o data/TempestRooted1327_WGSRate_OutRem.date.nwk -e 5 -s 10860 -f 1000 -w data/rate.txt
all.result=lsd2(inputTree="processed_data/wgs_alignment/all_rooted.nwk", inputDate ="processed_data/wgs_alignment/all.dates.txt", outFile="processed_data/wgs_alignment/all_wgsrate_lsd_CI", seqLen=11925, ZscoreOutlier=5, givenRate=rate, confidenceInterval=1000)
# rate 0.000203983, tMRCA 1894.77 , objective function 0.21766

