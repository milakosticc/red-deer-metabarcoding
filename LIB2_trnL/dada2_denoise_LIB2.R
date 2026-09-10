library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_trnLc-trnLh"

filtered <- file.path(base, "qualFiltered_out")
err_dir  <- file.path(base, "dada2_out/error_models")
out      <- file.path(base, "dada2_out/denoised")

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

sample.names.F <- sub("_filtered_R1.fastq.gz$", "", basename(fnFs))
sample.names.R <- sub("_filtered_R2.fastq.gz$", "", basename(fnRs))

stopifnot(
  length(fnFs) == length(fnRs),
  identical(sample.names.F, sample.names.R)
)

sample.names <- sample.names.F

errF <- readRDS(
  file.path(err_dir, "LIB2_error_model_R1.rds")
)

errR <- readRDS(
  file.path(err_dir, "LIB2_error_model_R2.rds")
)

cat("Denoising", length(sample.names), "samples\n")

cat("\nDenoising R1...\n")
dadaFs <- dada(
  fnFs,
  err = errF,
  multithread = TRUE
)

cat("\nDenoising R2...\n")
dadaRs <- dada(
  fnRs,
  err = errR,
  multithread = TRUE
)

names(dadaFs) <- sample.names
names(dadaRs) <- sample.names

saveRDS(
  dadaFs,
  file.path(out, "LIB2_dada_R1.rds")
)

saveRDS(
  dadaRs,
  file.path(out, "LIB2_dada_R2.rds")
)

getN <- function(x) sum(getUniques(x))

denoise_summary <- data.frame(
  Sample = sample.names,
  Denoised_R1 = sapply(dadaFs, getN),
  Denoised_R2 = sapply(dadaRs, getN)
)

write.csv(
  denoise_summary,
  file.path(out, "LIB2_denoising_summary.csv"),
  row.names = FALSE
)

print(denoise_summary)

cat("\nDenoising completed successfully.\n")
