library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

asv_dir <- file.path(base, "dada2_out/asv_table")
out <- file.path(base, "dada2_out/nonchim")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Load raw ASV table
seqtab <- readRDS(
  file.path(asv_dir, "LIB1_raw_ASV_table.rds")
)

asv_before <- ncol(seqtab)
reads_before <- sum(seqtab)

cat("ASVs before chimera removal:", asv_before, "\n")
cat("Abundance before chimera removal:", reads_before, "\n")

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

asv_after <- ncol(seqtab.nochim)
reads_after <- sum(seqtab.nochim)

asv_retained_pct <- 100 * asv_after / asv_before
reads_retained_pct <- 100 * reads_after / reads_before

cat("\nASVs after chimera removal:", asv_after, "\n")
cat("ASVs retained:", round(asv_retained_pct, 2), "%\n")

cat("\nAbundance after chimera removal:", reads_after, "\n")
cat("Abundance retained:", round(reads_retained_pct, 2), "%\n")

# Save non-chimeric ASV table
saveRDS(
  seqtab.nochim,
  file.path(out, "LIB1_nonchim_ASV_table.rds")
)

write.csv(
  seqtab.nochim,
  file.path(out, "LIB1_nonchim_ASV_table.csv")
)

# Save summary
summary <- data.frame(
  ASVs_before = asv_before,
  ASVs_after = asv_after,
  ASV_retained_percent = asv_retained_pct,
  Abundance_before = reads_before,
  Abundance_after = reads_after,
  Abundance_retained_percent = reads_retained_pct
)

write.csv(
  summary,
  file.path(out, "LIB1_chimera_summary.csv"),
  row.names = FALSE
)

cat("\nChimera removal completed successfully.\n")
