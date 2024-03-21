
#' merge_sequences_by_id function ----------------------------------------------------------------------
#' @author
#' Kirstyn Brunker
#' written on: 15 Dec 2023
#' @description
#' Will combine sequence records that have the same isolate id (e.g. multi-gene)
#' @details
#' Requires the library seqinr
#' @example
#' input_file <- "seq.fasta" ; merged_sequences <- merge_sequences_by_id(input_file)

merge_sequences_by_id <- function(file_path, fixids = FALSE) {
  # Read the FASTA file
  sequences <- readDNAStringSet(file_path)
  
  # Fix IDs if fixids is set to TRUE
  if (fixids) {
    names(sequences) <- gsub("/[0-9]+-[0-9]+$", "", names(sequences))
  }
  
  # Identify unique sequence IDs
  unique_ids <- unique(names(sequences))
  merged_sequences <- DNAStringSet(rep("", length(unique_ids)))
  
  # Merge sequences by column and replace gaps ("-") with letters if present
  for (id in unique_ids) {
    same_id_seqs <- as.matrix(sequences[names(sequences) == id])
    merged_seq <- apply(same_id_seqs, 2, function(col) paste0(col[col != "-"], collapse = ""))
    merged_seq2 <- gsub("^$", "-", merged_seq)
    merged_seq2 <- paste(merged_seq2, collapse = "")
    merged_sequences[unique_ids == id] <- DNAStringSet(paste(merged_seq2, collapse = ""), use.names = TRUE)
  }
  names(merged_sequences) <- unique_ids
  return(merged_sequences)
}


