library(data.table)

# Load chimera-filtered ASV table
seqtab <- readRDS("dada2_out/seqtab_nochim.rds")

# Bioscanflow UNCROSS2 parameters
set_f <- 0.03
set_p <- 1
set_tmin <- 0.1

# UNCROSS score function from Bioscanflow
uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1) {

  z <- f * N / n
  sc <- 2 / (1 + exp(x / z)^p)

  res <- data.table(
    Score = sc,
    TagJump = sc >= tmin
  )

  return(res)
}

# --------------------------------
# Convert ASV table to long format
# --------------------------------

dt <- as.data.table(
  seqtab,
  keep.rownames = "Sample"
)

long <- melt(
  dt,
  id.vars = "Sample",
  variable.name = "ASV_sequence",
  value.name = "Abundance"
)

# Only observed ASV/sample combinations
long <- long[Abundance > 0]

# Total abundance of each ASV across the LIB3 pool
long[, Total_ASV_abundance := sum(Abundance), by = ASV_sequence]

# Number of samples in the pool
n_samples <- nrow(seqtab)

# Calculate UNCROSS2 score
scores <- uncross_score(
  x = long$Abundance,
  N = long$Total_ASV_abundance,
  n = n_samples,
  f = set_f,
  tmin = set_tmin,
  p = set_p
)

long[, Score := scores$Score]
long[, TagJump := scores$TagJump]

# Save complete UNCROSS2 assessment
fwrite(
  long,
  "dada2_out/UNCROSS2_scores_LIB3.csv"
)

# --------------------------------
# Remove tag-jump abundances
# --------------------------------

long[TagJump == TRUE, Abundance := 0]

# Convert back to ASV abundance table
filtered.long <- long[, .(Abundance = sum(Abundance)),
                      by = .(Sample, ASV_sequence)]

filtered.wide <- dcast(
  filtered.long,
  Sample ~ ASV_sequence,
  value.var = "Abundance",
  fill = 0
)

sample.names <- filtered.wide$Sample
filtered.wide[, Sample := NULL]

seqtab.uncross <- as.matrix(filtered.wide)
rownames(seqtab.uncross) <- sample.names

# Preserve original sample order
seqtab.uncross <- seqtab.uncross[rownames(seqtab), , drop = FALSE]

# Remove ASVs that became zero in every sample
seqtab.uncross <- seqtab.uncross[
  ,
  colSums(seqtab.uncross) > 0,
  drop = FALSE
]

# Save result
saveRDS(
  seqtab.uncross,
  "dada2_out/seqtab_uncross2.rds"
)

write.csv(
  seqtab.uncross,
  "dada2_out/ASV_table_uncross2_sequences.csv",
  quote = FALSE
)

# --------------------------------
# Summary
# --------------------------------

cat("\n=== UNCROSS2 FILTERING ===\n")

cat("Samples:", nrow(seqtab.uncross), "\n")

cat("ASVs before:", ncol(seqtab), "\n")
cat("ASVs after:", ncol(seqtab.uncross), "\n")

cat("Reads before:", sum(seqtab), "\n")
cat("Reads after:", sum(seqtab.uncross), "\n")

cat(
  "ASV retention:",
  round(100 * ncol(seqtab.uncross) / ncol(seqtab), 2),
  "%\n"
)

cat(
  "Read abundance retention:",
  round(100 * sum(seqtab.uncross) / sum(seqtab), 2),
  "%\n"
)

cat(
  "Tag-jump sample-ASV occurrences removed:",
  sum(long$TagJump),
  "\n"
)
