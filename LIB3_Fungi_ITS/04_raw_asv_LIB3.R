library(dada2)

# Load merged reads
mergers <- readRDS("dada2_out/mergers.rds")

# Construct raw ASV sequence table
seqtab <- makeSequenceTable(mergers)

# Save intermediate output
saveRDS(
  seqtab,
  "dada2_out/seqtab_raw.rds"
)

write.csv(
  seqtab,
  "dada2_out/ASV_table_raw_sequences.csv",
  quote = FALSE
)

# Summary
asv.lengths <- nchar(colnames(seqtab))

cat("\n=== RAW ASV TABLE ===\n")
cat("Samples:", nrow(seqtab), "\n")
cat("ASVs:", ncol(seqtab), "\n")
cat("Total abundance:", sum(seqtab), "\n")
cat("Minimum ASV length:", min(asv.lengths), "bp\n")
cat("Maximum ASV length:", max(asv.lengths), "bp\n")

cat("\nASV length distribution:\n")
print(sort(table(asv.lengths), decreasing = TRUE))
