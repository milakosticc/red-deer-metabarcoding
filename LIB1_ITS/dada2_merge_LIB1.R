library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

filtered <- file.path(base, "qualFiltered_out")
denoised <- file.path(base, "dada2_out/denoised")
out      <- file.path(base, "dada2_out/merged")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Input filtered reads
fnFs <- sort(list.files(
  filtered,
  pattern = "_filtered_R1.fastq.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  filtered,
  pattern = "_filtered_R2.fastq.gz$",
  full.names = TRUE
))

sample.names <- sub("_filtered_R1.fastq.gz$", "", basename(fnFs))

# Load denoised DADA2 objects
dadaFs <- readRDS(
  file.path(denoised, "LIB1_dada_R1.rds")
)

dadaRs <- readRDS(
  file.path(denoised, "LIB1_dada_R2.rds")
)

# Merge paired reads
mergers <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  verbose = TRUE
)

names(mergers) <- sample.names

# Save intermediate merged object
saveRDS(
  mergers,
  file.path(out, "LIB1_mergers.rds")
)

# Count successfully merged reads
getN <- function(x) sum(x$abundance)

merged_summary <- data.frame(
  Sample = sample.names,
  Merged_reads = sapply(mergers, getN)
)

write.csv(
  merged_summary,
  file.path(out, "LIB1_merged_summary.csv"),
  row.names = FALSE
)

print(merged_summary)

cat(
  "\nTotal successfully merged read pairs:",
  sum(merged_summary$Merged_reads),
  "\n"
)

cat("\nPaired-end merging completed successfully.\n")
