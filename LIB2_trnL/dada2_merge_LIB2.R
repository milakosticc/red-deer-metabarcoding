library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_trnLc-trnLh"

filtered <- file.path(base, "qualFiltered_out")
denoised <- file.path(base, "dada2_out/denoised")
out      <- file.path(base, "dada2_out/merged")

dir.create(out, recursive = TRUE, showWarnings = FALSE)

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

dadaFs <- readRDS(
  file.path(denoised, "LIB2_dada_R1.rds")
)

dadaRs <- readRDS(
  file.path(denoised, "LIB2_dada_R2.rds")
)

mergers <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  verbose = TRUE
)

names(mergers) <- sample.names

saveRDS(
  mergers,
  file.path(out, "LIB2_mergers.rds")
)

getN <- function(x) sum(x$abundance)

merged_summary <- data.frame(
  Sample = sample.names,
  Merged_reads = sapply(mergers, getN)
)

write.csv(
  merged_summary,
  file.path(out, "LIB2_merged_summary.csv"),
  row.names = FALSE
)

print(merged_summary)

cat(
  "\nTotal successfully merged read pairs:",
  sum(merged_summary$Merged_reads),
  "\n"
)

cat("\nPaired-end merging completed successfully.\n")
