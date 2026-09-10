base <- "/home/mila/aapraksa/samples/working/Plants_trnLc-trnLh"

# Previous preprocessing summary
pre <- read.csv(
  file.path(base, "seq_count_summary.csv"),
  stringsAsFactors = FALSE
)

# Denoising
den <- read.csv(
  file.path(base, "dada2_out/denoised/LIB2_denoising_summary.csv"),
  stringsAsFactors = FALSE
)

# Merging
mer <- read.csv(
  file.path(base, "dada2_out/merged/LIB2_merged_summary.csv"),
  stringsAsFactors = FALSE
)

# Non-chimeric table
nochim <- readRDS(
  file.path(base, "dada2_out/nonchim/LIB2_nonchim_ASV_table.rds")
)

# UNCROSS2 table
uncross <- readRDS(
  file.path(base, "dada2_out/uncross2/LIB2_UNCROSS2_filtered_ASV_table.rds")
)

# Totals
input_reads <- sum(pre$reads_in)
filtered_reads <- sum(pre$reads_out)

den_R1 <- sum(den$Denoised_R1)
den_R2 <- sum(den$Denoised_R2)

merged_reads <- sum(mer$Merged_reads)

nonchim_reads <- sum(nochim)
uncross_reads <- sum(uncross)

tracking <- data.frame(
  Step = c(
    "Input",
    "Quality-filtered",
    "Denoised R1",
    "Denoised R2",
    "Merged",
    "Non-chimeric",
    "UNCROSS2-filtered",
    "Contaminant-filtered"
  ),

  Reads = c(
    input_reads,
    filtered_reads,
    den_R1,
    den_R2,
    merged_reads,
    nonchim_reads,
    uncross_reads,
    NA
  ),

  Status = c(
    "Completed",
    "Completed",
    "Completed",
    "Completed",
    "Completed",
    "Completed",
    "Completed",
    "Not performed - no suitable controls/metadata"
  )
)

write.csv(
  tracking,
  file.path(base, "dada2_out/final/LIB2_final_sequence_tracking.csv"),
  row.names = FALSE
)

# ASV tracking
asv_tracking <- data.frame(
  Step = c(
    "Raw ASV table",
    "Non-chimeric",
    "UNCROSS2-filtered",
    "Final"
  ),

  ASVs = c(
    405,
    ncol(nochim),
    ncol(uncross),
    ncol(uncross)
  )
)

write.csv(
  asv_tracking,
  file.path(base, "dada2_out/final/LIB2_final_ASV_tracking.csv"),
  row.names = FALSE
)

cat("=== FINAL LIB2 SEQUENCE TRACKING ===\n")
print(tracking)

cat("\n=== FINAL LIB2 ASV TRACKING ===\n")
print(asv_tracking)

cat("\nFinal tracking completed successfully.\n")
