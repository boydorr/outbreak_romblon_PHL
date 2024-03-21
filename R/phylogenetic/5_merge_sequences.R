
#'---------------------------------------------------------
#'title: Merging Sequences 
#'date: 19/12/2023
#'---------------------------------------------------------
# Scripts to match fasta sequence names to the metadata isolate ID for cleaning data
rm(list=ls())  # This clears anything currently in the environment

#############################################
#  INSTALL PACKAGES                         #
############################################# 
library(seqinr)
library(Biostrings)

#############################################
#              LOAD FUNCTIONS               #
#############################################
source("scripts/prepareSequences_functions.R")
# These scripts were created by Kirstyn Brunker 15/Dec/2023

# edit below for specific input and output files

# assign the input file
input_file <- "processed_data/concatenated_alignment/ph_all_combined_690.aln.fasta"
# apply the function
merged_sequences <- merge_sequences_by_id(input_file, fixids=T)
nseq=length(merged_sequences)
# save the output
writeXStringSet(merged_sequences, paste0("processed_data/concatenated_alignment/ph_all_merged_by_id_",nseq,".fasta"), format = "fasta")

# and generate associated metadata file
meta=read.csv("processed_data/data_prep/metadata/ph_metadata_693_clean.csv")

# keep the record with the most non-NA values across all columns
collapsed_meta <- meta %>%
  group_by(Sample.ID) %>%
  summarise_all(~ if_else(all(is.na(.)), NA, .[which.max(!is.na(.))])) %>%
  ungroup()
# save the output
write.csv(collapsed_meta , paste0("processed_data/concatenated_alignment/ph_metadata_merged_by_id_",nseq,".csv"), row.names=F)
collapsed_meta$Sample.ID[!collapsed_meta$Sample.ID %in% names(merged_sequences)]

# output in format for iqtree lsd
lsd=cbind(meta$Sample.ID, meta$decimal_date)
write.table(lsd,paste0("processed_data/concatenated_alignment/ph_metadata_merged_by_id_",nseq,"_lsd_dates.txt"), quote=F, sep="\t", col.names=F, row.names=F)

# and for pastml analysis
pastml=data.frame(Sample.ID=meta$Sample.ID, Region=meta$Region, Province=meta$Province)
write.csv(pastml,paste0("processed_data/concatenated_alignment/ph_metadata_merged_by_id_",nseq,"_pastml.csv"),row.names=F)
       