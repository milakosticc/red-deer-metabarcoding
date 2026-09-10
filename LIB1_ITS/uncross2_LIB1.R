library(data.table)

base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

input_file <- file.path(
  base,
  "dada2_out/nonchim/LIB1_nonchim_ASV_table.rds"
)

out <- file.path(
  base,
  "dada2_out/uncross2"
)

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# UNCROSS2 parameters used in BioScanFlow
set_f <- 0.03
set_p <- 1

# Load non-chimeric ASV table
tab <- readRDS(input_file)

cat("=== UNCROSS2 INPUT ===\n")
cat("Samples:", nrow(tab), "\n")
cat("ASVs:", ncol(tab), "\n")
cat("Total abundance:", sum(tab), "\n")

asvs_before <- ncol(tab)
reads_before <- sum(tab)

# Convert to long format
ASVTABW <- as.data.table(
  t(tab),
  keep.rownames = TRUE
)

colnames(ASVTABW)[1] <- "ASV"

ASVTAB <- melt(
  ASVTABW,
  id.vars = "ASV",
  variable.name = "SampleID",
  value.name = "Abundance"
)

# Keep only observed ASV/sample combinations
ASVTAB <- ASVTAB[Abundance > 0]

# Total abundance of each ASV across LIB1
ASVTAB[
  ,
  Total := sum(Abundance, na.rm = TRUE),
  by = ASV
]

# UNCROSS2 score function 
uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1) {

  z <- f * N / n

  sc <- 2 / (1 + exp(x / z)^p)

  res <- data.table(
    Score = sc,
    TagJump = sc >= tmin
  )

  return(res)
}

n_samples <- length(unique(ASVTAB$SampleID))

scores <- uncross_score(
  x = ASVTAB$Abundance,
  N = ASVTAB$Total,
  n = n_samples,
  f = set_f,
  p = set_p
)

ASVTAB <- cbind(
  ASVTAB,
  scores
)

# Save diagnostic score table
fwrite(
  ASVTAB,
  file.path(out, "LIB1_UNCROSS2_scores.tsv"),
  sep = "\t"
)

# Tag-jump statistics
tagjump_events <- sum(
  ASVTAB$TagJump,
  na.rm = TRUE
)

tagjump_reads <- sum(
  ASVTAB[TagJump == TRUE]$Abundance,
  na.rm = TRUE
)

cat("\n=== UNCROSS2 TAG-JUMP RESULTS ===\n")
cat("Tag-jump events:", tagjump_events, "\n")
cat("Tag-jump reads:", tagjump_reads, "\n")
cat(
  "Read abundance removed:",
  round(100 * tagjump_reads / reads_before, 2),
  "%\n"
)

# Remove tag-jump observations
ASVTAB.filtered <- ASVTAB[
  TagJump == FALSE
]

# Convert back to ASV abundance table
RES <- dcast(
  ASVTAB.filtered,
  ASV ~ SampleID,
  value.var = "Abundance",
  fill = 0
)

asv_sequences <- RES$ASV

filtered_matrix <- as.matrix(
  RES[, -"ASV"]
)

rownames(filtered_matrix) <- asv_sequences

seqtab.uncross <- t(filtered_matrix)

storage.mode(seqtab.uncross) <- "numeric"

# Keep original sample order
seqtab.uncross <- seqtab.uncross[
  rownames(tab),
  ,
  drop = FALSE
]

asvs_after <- ncol(seqtab.uncross)
reads_after <- sum(seqtab.uncross)

# Save filtered ASV table
saveRDS(
  seqtab.uncross,
  file.path(out, "LIB1_UNCROSS2_filtered_ASV_table.rds")
)

write.csv(
  seqtab.uncross,
  file.path(out, "LIB1_UNCROSS2_filtered_ASV_table.csv")
)

# Summary
summary <- data.frame(
  ASVs_before = asvs_before,
  ASVs_after = asvs_after,
  ASV_retained_percent =
    100 * asvs_after / asvs_before,

  Abundance_before = reads_before,
  Abundance_after = reads_after,
  Abundance_retained_percent =
    100 * reads_after / reads_before,

  TagJump_events = tagjump_events,
  TagJump_reads_removed = tagjump_reads,
  TagJump_read_percent_removed =
    100 * tagjump_reads / reads_before,

  f_parameter = set_f,
  p_parameter = set_p
)

write.csv(
  summary,
  file.path(out, "LIB1_UNCROSS2_summary.csv"),
  row.names = FALSE
)

cat("\n=== UNCROSS2 OUTPUT ===\n")
cat("ASVs before:", asvs_before, "\n")
cat("ASVs after:", asvs_after, "\n")
cat(
  "ASVs retained:",
  round(100 * asvs_after / asvs_before, 2),
  "%\n"
)

cat("\nAbundance before:", reads_before, "\n")
cat("Abundance after:", reads_after, "\n")
cat(
  "Abundance retained:",
  round(100 * reads_after / reads_before, 2),
  "%\n"
)

cat("\nUNCROSS2 filtering completed successfully.\n")
