DNAmotifs <- function(fasta_file, ws, cut_off) {

  # Function to read a multi-FASTA file
  read_multi_fasta <- function(file_path) {
    dna_sequences <- readDNAStringSet(file_path)
    return(dna_sequences)
  }

  # Read the multi-FASTA file
  sequences <- read_multi_fasta(fasta_file)

  # Initialize a list to store results
  results <- list()

  # Initialize a vector to store motifs
  all_motifs <- character(0)

  # Loop through each sequence
  for (i in seq_along(sequences)) {
    # Extract sequence and sequence ID
    seq <- as.character(sequences[[i]])
    seq_id <- names(sequences)[i]

    # Call the C++ function for motif analysis
    result <- .generate_frequency_table(seq, seq_id, nchar(seq), ws, cut_off,
                                       all_motifs, length(all_motifs))

    # Extract motifs and counts from the result
    motifs <- as.character(result$motif)
    counts <- as.integer(result$count)

    # Update the all_motifs vector
    all_motifs <- unique(c(all_motifs, motifs))

    # Store the counts in a named vector
    counts_named <- setNames(counts, motifs)

    # Append the named vector to results
    results[[seq_id]] <- counts_named
  }

  # Create a final data frame with sequence IDs in the first column and motifs
  # as column names
  motif_columns <- unique(unlist(lapply(results, names)))
  final_results <- data.frame(matrix(0, nrow = length(sequences),
                                     ncol = length(motif_columns) + 1))
  colnames(final_results) <- c("sequence_id", motif_columns)

  # Populate the final data frame
  for (i in seq_along(sequences)) {
    seq_id <- names(sequences)[i]
    final_results[i, "sequence_id"] <- seq_id
    if (!is.null(results[[seq_id]])) {
      for (motif in names(results[[seq_id]])) {
        final_results[i, motif] <- results[[seq_id]][motif]
      }
    }
  }

  # Return the final results
  return(final_results)
}
