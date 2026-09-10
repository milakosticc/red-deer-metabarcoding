base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

input_file <- file.path(
  base,
  "dada2_out/uncross2/LIB1_UNCROSS2_filtered_ASV_table.rds"
)

out <- file.path(base, "dada2_out/final")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load final cleaned ASV table
tab <- readRDS(input_file)

# Original column names are nucleotide sequences
sequences <- colnames(tab)

# Assign ASV IDs
asv_ids <- paste0("ASV", seq_along(sequences))

final_tab <- tab
colnames(final_tab) <- asv_ids

# Save final abundance table
saveRDS(
  final_tab,
  file.path(out, "LIB1_final_ASV_table.rds")
)

write.csv(
  final_tab,
  file.path(out, "LIB1_final_ASV_table.csv")
)

# ASV ID - sequence mapping
mapping <- data.frame(
  ASV_ID = asv_ids,
  Sequence = sequences,
  Total_abundance = colSums(tab),
  stringsAsFactors = FALSE
)

write.csv(
  mapping,
  file.path(out, "LIB1_final_ASV_mapping.csv"),
  row.names = FALSE
)

# FASTA export
fasta <- as.vector(
  rbind(
    paste0(">", asv_ids),
    sequences
  )
)

writeLines(
  fasta,
  file.path(out, "LIB1_final_ASVs.fasta")
)

# Verify IDs
stopifnot(
  identical(colnames(final_tab), mapping$ASV_ID),
  length(sequences) == ncol(final_tab)
)

cat("=== FINAL LIB1 ASV OUTPUT ===\n")
cat("Samples:", nrow(final_tab), "\n")
cat("ASVs:", ncol(final_tab), "\n")
cat("Total abundance:", sum(final_tab), "\n")

cat("\nFirst ASV IDs:\n")
print(head(colnames(final_tab), 10))

cat("\nFASTA and abundance-table ASV identifiers correspond correctly.\n")
cat("Final LIB1 ASV export completed successfully.\n")
