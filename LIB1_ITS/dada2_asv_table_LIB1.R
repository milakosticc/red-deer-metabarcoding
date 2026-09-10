library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

merged_dir <- file.path(base, "dada2_out/merged")
out <- file.path(base, "dada2_out/asv_table")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load merged paired-end reads
mergers <- readRDS(
  file.path(merged_dir, "LIB1_mergers.rds")
)

# Construct raw ASV sequence table
seqtab <- makeSequenceTable(mergers)

cat("Number of samples:", nrow(seqtab), "\n")
cat("Number of ASVs:", ncol(seqtab), "\n")
cat("Total abundance:", sum(seqtab), "\n")

# Save raw ASV table
saveRDS(
  seqtab,
  file.path(out, "LIB1_raw_ASV_table.rds")
)

write.csv(
  seqtab,
  file.path(out, "LIB1_raw_ASV_table.csv")
)

# ASV length distribution
asv_lengths <- nchar(colnames(seqtab))
length_summary <- table(asv_lengths)

write.csv(
  as.data.frame(length_summary),
  file.path(out, "LIB1_ASV_length_distribution.csv"),
  row.names = FALSE
)

cat("\nASV length distribution:\n")
print(length_summary)

cat("\nRaw ASV table created successfully.\n")
