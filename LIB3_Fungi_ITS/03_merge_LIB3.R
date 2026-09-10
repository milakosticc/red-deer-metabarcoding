library(dada2)

# Load intermediate DADA2 objects
derepFs <- readRDS("dada2_out/derepFs.rds")
derepRs <- readRDS("dada2_out/derepRs.rds")

dadaFs <- readRDS("dada2_out/dadaFs.rds")
dadaRs <- readRDS("dada2_out/dadaRs.rds")

# Merge paired-end reads
mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = TRUE
)

# Keep intermediate output
saveRDS(
  mergers,
  "dada2_out/mergers.rds"
)

getN <- function(x) sum(getUniques(x))

denoisedF <- sapply(dadaFs, getN)
denoisedR <- sapply(dadaRs, getN)
merged <- sapply(mergers, getN)

# Number of possible paired reads per sample
candidatePairs <- pmin(denoisedF, denoisedR)

merge.track <- data.frame(
  sample = names(merged),
  denoisedF = denoisedF,
  denoisedR = denoisedR,
  candidatePairs = candidatePairs,
  merged = merged,
  merged_percent = round(100 * merged / candidatePairs, 2)
)

write.csv(
  merge.track,
  "dada2_out/merging_counts.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== MERGING COUNTS ===\n")
print(merge.track)

cat("\nTOTAL CANDIDATE PAIRS:", sum(candidatePairs), "\n")
cat("TOTAL MERGED:", sum(merged), "\n")
cat(
  "OVERALL MERGING RETENTION:",
  round(100 * sum(merged) / sum(candidatePairs), 2),
  "%\n"
)
