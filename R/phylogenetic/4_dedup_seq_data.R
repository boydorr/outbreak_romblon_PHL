### deduplicate the glue sequences
## i.e. remove partial seq if a WGS exists or if a newer version exists
library(seqinr)
#library(sequences)

glue=read.fasta("processed_data/data_prep/sequences/ph_rabv_glue_wg_isolateIds.fasta")
wgs=read.fasta("processed_data/data_prep/sequences/ph_concat_wgs_78.aln.fasta")
glue.meta=read.csv("raw_data/rabv-glue_data/ph_metadata_615.csv")

# identify duplicates
duplicate.names=names(glue)[match(names(wgs),names(glue),nomatch = 0)]

## remove from the concatenated fasta
glue.dedup=glue[-match(names(wgs),names(glue),nomatch = 0)]
## and remove from the metadata
glue.meta2 = glue.meta[!(glue.meta$sequence.isolate %in% duplicate.names), ]

## remove duplicate isolate sequences that have a WGS representative- WGS trumps and replaces any partial seq
# Find duplicate IDs
duplicate <- glue.meta[duplicated(glue.meta$sequence.isolate) | duplicated(glue.meta$sequence.isolate, fromLast = TRUE),]

## assess each group of duplicates and remove any duplicates if there is at least one WGS representative
filter_group <- function(duplicate) {
  # Check if any record has sequence length > 11700
  has_large_sequence <- any(duplicate$sequence.gb_length > 11700)
  
  # If yes, keep only the one with sequence length > 11700, otherwise keep all
  if (has_large_sequence) {
    return(duplicate[duplicate$sequence.gb_length > 11700, ])
    print(duplicate$sequence.sequenceID)
  } else {
    return(duplicate)
  }
}

# Apply the filtering function to each group of duplicates
filtered_data <- do.call(rbind, lapply(split(glue.meta2, glue.meta2$sequence.isolate), filter_group))
write.csv(filtered_data, "processed_data/data_prep/metadata/ph_metadata_612.csv")
# names of duplicated sequence isolates removed
removed <- unique(glue.meta2$sequence.isolate[!glue.meta2$sequence.sequenceID %in% filtered_data$sequence.sequenceID])


# remove from fasta file
# Initialize a vector to store the results
sequence_lengths <- numeric(length(glue.dedup))

# Loop through each sequence in the multifasta file
for (i in seq_along(glue.dedup)) {
  # Remove dashes and calculate sequence length
  sequence_lengths[i] <- sum(nchar(gsub("-", "", glue.dedup[i][[1]])))
}

# Save the original index before removal
original_index <- which(names(glue.dedup) == removed)
# Perform the removal operation
identify <- which(names(glue.dedup) == removed)
to_remove <- which(sequence_lengths[identify] < 11700)

# Adjust the original index based on the removal operation
to_remove_original_index <- original_index[to_remove]
# remove from fasta
glue.dedup2=glue.dedup[-to_remove_original_index]
write.fasta(glue.dedup2, names=names(glue.dedup2),"processed_data/data_prep/sequences/ph_rabv_glue_wg_isolateIds_dedup.fasta")            
