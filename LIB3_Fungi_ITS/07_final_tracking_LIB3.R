library(dada2)

# -----------------------------
# INPUT + QUALITY FILTERING
# -----------------------------

filter_counts <- read.csv("seq_count_summary.csv")

input_reads <- sum(filter_counts$reads_in)
filtered_reads <- sum(filter_counts$reads_out)

# -----------------------------
# DENOISING
# -----------------------------

denoising <- read.csv("dada2_out/denoising_counts.csv")

denoised_R1 <- sum(denoising$denoisedF)
denoised_R2 <- sum(denoising$denoisedR)

# -----------------------------
# MERGING
# -----------------------------

merging <- read.csv("dada2_out/merging_counts.csv")
merged_reads <- sum(merging$merged)

# -----------------------------
# ASV TABLES
# -----------------------------

seqtab.raw <- readRDS("dada2_out/seqtab_raw.rds")
seqtab.nochim <- readRDS("dada2_out/seqtab_nochim.rds")
seqtab.uncross <- readRDS("dada2_out/seqtab_uncross2.rds")

nonchim_reads <- sum(seqtab.nochim)
uncross_reads <- sum(seqtab.uncross)

# -----------------------------
# SEQUENCE TRACKING
# -----------------------------

sequence_tracking <- data.frame(
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
    denoised_R1,
    denoised_R2,
    merged_reads,
    nonchim_reads,
    uncross_reads,
    NA
  ),
  Status = c(
    rep("Completed", 7),
    "Not performed - no suitable controls/metadata"
  )
)

# -----------------------------
# ASV TRACKING
# -----------------------------

asv_tracking <- data.frame(
  Step = c(
    "Raw ASV table",
    "Non-chimeric",
    "UNCROSS2-filtered",
    "Final"
  ),
  ASVs = c(
    ncol(seqtab.raw),
    ncol(seqtab.nochim),
    ncol(seqtab.uncross),
    ncol(seqtab.uncross)
  )
)

# -----------------------------
# SAVE
# -----------------------------

dir.create(
  "dada2_out/final",
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  sequence_tracking,
  "dada2_out/final/LIB3_final_sequence_tracking.csv",
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  asv_tracking,
  "dada2_out/final/LIB3_final_ASV_tracking.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("=== FINAL LIB3 SEQUENCE TRACKING ===\n")
print(sequence_tracking)

cat("\n=== FINAL LIB3 ASV TRACKING ===\n")
print(asv_tracking)

cat("\nFinal tracking completed successfully.\n")
