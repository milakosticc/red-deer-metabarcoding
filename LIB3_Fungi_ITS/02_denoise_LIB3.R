library(dada2)

path <- "qualFiltered_out"

fnFs <- sort(list.files(
  path,
  pattern = "_filtered_R1.fastq.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  path,
  pattern = "_filtered_R2.fastq.gz$",
  full.names = TRUE
))

sample.names <- sub("_filtered_R1.fastq.gz$", "", basename(fnFs))
sample.names.R <- sub("_filtered_R2.fastq.gz$", "", basename(fnRs))

stopifnot(length(fnFs) == length(fnRs))
stopifnot(all(sample.names == sample.names.R))

names(fnFs) <- sample.names
names(fnRs) <- sample.names

# Load learned error models
errF <- readRDS("dada2_out/error_model_F.rds")
errR <- readRDS("dada2_out/error_model_R.rds")

# Dereplicate
derepFs <- derepFastq(fnFs)
derepRs <- derepFastq(fnRs)

names(derepFs) <- sample.names
names(derepRs) <- sample.names

saveRDS(derepFs, "dada2_out/derepFs.rds")
saveRDS(derepRs, "dada2_out/derepRs.rds")

# Denoise
dadaFs <- dada(
  derepFs,
  err = errF,
  multithread = TRUE
)

dadaRs <- dada(
  derepRs,
  err = errR,
  multithread = TRUE
)

saveRDS(dadaFs, "dada2_out/dadaFs.rds")
saveRDS(dadaRs, "dada2_out/dadaRs.rds")

# Count reads after denoising
getN <- function(x) sum(getUniques(x))

denoised <- data.frame(
  sample = sample.names,
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN)
)

write.csv(
  denoised,
  "dada2_out/denoising_counts.csv",
  row.names = FALSE,
  quote = FALSE
)

cat("\n=== DENOISING COUNTS ===\n")
print(denoised)

cat("\nTOTAL R1 DENOISED:", sum(denoised$denoisedF), "\n")
cat("TOTAL R2 DENOISED:", sum(denoised$denoisedR), "\n")
cat(
  "R1 RETAINED:",
  round(100 * sum(denoised$denoisedF) / 351497, 2),
  "%\n"
)
cat(
  "R2 RETAINED:",
  round(100 * sum(denoised$denoisedR) / 351497, 2),
  "%\n"
)
