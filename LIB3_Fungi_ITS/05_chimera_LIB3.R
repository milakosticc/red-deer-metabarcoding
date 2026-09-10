library(dada2)

# Load raw ASV table
seqtab <- readRDS("dada2_out/seqtab_raw.rds")

# Remove chimeric ASVs
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Save chimera-filtered table
saveRDS(
  seqtab.nochim,
  "dada2_out/seqtab_nochim.rds"
)

write.csv(
  seqtab.nochim,
  "dada2_out/ASV_table_nochim_sequences.csv",
  quote = FALSE
)

# Summary
asv.before <- ncol(seqtab)
asv.after <- ncol(seqtab.nochim)

reads.before <- sum(seqtab)
reads.after <- sum(seqtab.nochim)

cat("\n=== CHIMERA REMOVAL ===\n")

cat("ASVs before:", asv.before, "\n")
cat("ASVs after:", asv.after, "\n")
cat("Chimeric ASVs removed:", asv.before - asv.after, "\n")

cat(
  "ASV retention:",
  round(100 * asv.after / asv.before, 2),
  "%\n"
)

cat("\nReads before:", reads.before, "\n")
cat("Reads after:", reads.after, "\n")

cat(
  "Read abundance retention:",
  round(100 * reads.after / reads.before, 2),
  "%\n"
)
