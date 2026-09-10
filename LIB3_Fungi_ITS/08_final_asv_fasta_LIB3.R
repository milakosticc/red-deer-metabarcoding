# Load final UNCROSS2-filtered table
seqtab <- readRDS("dada2_out/seqtab_uncross2.rds")

# Original DNA sequences are currently column names
asv_sequences <- colnames(seqtab)

# Create simple ASV identifiers
asv_ids <- paste0("ASV", seq_along(asv_sequences))

# Rename columns in abundance table
colnames(seqtab) <- asv_ids

# Create output directory
dir.create(
  "dada2_out/final",
  showWarnings = FALSE,
  recursive = TRUE
)

# -----------------------------
# FINAL ASV ABUNDANCE TABLE
# -----------------------------

write.csv(
  seqtab,
  "dada2_out/final/LIB3_final_ASV_table.csv",
  quote = FALSE
)

# -----------------------------
# ASV ID ↔ SEQUENCE MAPPING
# -----------------------------

mapping <- data.frame(
  ASV = asv_ids,
  Sequence = asv_sequences
)

write.csv(
  mapping,
  "dada2_out/final/LIB3_final_ASV_mapping.csv",
  row.names = FALSE,
  quote = FALSE
)

# -----------------------------
# FASTA
# -----------------------------

fasta_lines <- as.vector(
  rbind(
    paste0(">", asv_ids),
    asv_sequences
  )
)

writeLines(
  fasta_lines,
  "dada2_out/final/LIB3_ASVs.fasta"
)

# -----------------------------
# VERIFY IDENTIFIERS
# -----------------------------

table_ids <- colnames(seqtab)
fasta_ids <- asv_ids

identical_ids <- identical(table_ids, fasta_ids)

cat("=== FINAL LIB3 ASV OUTPUT ===\n")
cat("Samples:", nrow(seqtab), "\n")
cat("ASVs:", ncol(seqtab), "\n")
cat("Total abundance:", sum(seqtab), "\n")
cat("FASTA sequences:", length(asv_sequences), "\n")
cat("FASTA IDs correspond to abundance-table IDs:", identical_ids, "\n")

if (!identical_ids) {
  stop("ERROR: FASTA and abundance-table identifiers do not correspond.")
}

cat("\nFinal ASV table, FASTA and mapping created successfully.\n")
