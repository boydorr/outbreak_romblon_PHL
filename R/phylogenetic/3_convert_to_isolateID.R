#'---------------------------------------------------------
#'title: Data Cleaning 
#'date: 8/12/2023
#'---------------------------------------------------------
# Scripts to match fasta sequence names to the metadata isolate ID for cleaning data
rm(list=ls())  # This clears anything currently in the environment

#############################################
#  INSTALL PACKAGES                         #
############################################# 
library(seqinr)

#############################################
#  IMPORT DATA                              #
#############################################

## import rabv-glue sequences

seq=read.fasta("raw_data/rabv-glue_data/ph_rabv_glue_wg_alignment.fasta")

# import metadata
metadata=read.csv("raw_data/rabv-glue_data/ph_metadata_615.csv")

# match seq names to metadata
match=match(names(seq), metadata$sequence.sequenceID)

# check for any not matched
names(seq)[!names(seq) %in% metadata$sequence.sequenceID]

#create names for beast 
new_ids=metadata$sequence.isolate[match]

# write new file with sequence names as isolate ids
write.fasta(seq, names=new_ids,"processed_data/data_prep/sequences/ph_rabv_glue_wg_isolateIds.fasta")

