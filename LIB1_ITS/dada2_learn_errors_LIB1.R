library(dada2)

# Input directory
path <- "/home/mila/aapraksa/samples/working/Plants_ITS/qualFiltered_out"

# Output directory
out <- "/home/mila/aapraksa/samples/working/Plants_ITS/dada2_out/error_models"

# Get forward and reverse filtered reads
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

# Extract sample names
sample.names.F <- sub("_filtered_R1.fastq.gz$", "", basename(fnFs))
sample.names.R <- sub("_filtered_R2.fastq.gz$", "", basename(fnRs))

# Check that R1 and R2 correspond correctly
stopifnot(
  length(fnFs) == length(fnRs),
  identical(sample.names.F, sample.names.R)
)

cat("Number of samples:", length(fnFs), "\n")
cat("Samples:\n")
print(sample.names.F)

# Learn forward-read error rates
cat("\nLearning R1 error rates...\n")
errF <- learnErrors(
  fnFs,
  multithread = TRUE,
  randomize = TRUE
)

# Learn reverse-read error rates
cat("\nLearning R2 error rates...\n")
errR <- learnErrors(
  fnRs,
  multithread = TRUE,
  randomize = TRUE
)

# Save error models
saveRDS(
  errF,
  file.path(out, "LIB1_error_model_R1.rds")
)

saveRDS(
  errR,
  file.path(out, "LIB1_error_model_R2.rds")
)

# Save error plots
pdf(
  file.path(out, "LIB1_error_rates_R1.pdf"),
  width = 8,
  height = 7
)
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf(
  file.path(out, "LIB1_error_rates_R2.pdf"),
  width = 8,
  height = 7
)
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("\nError learning completed successfully.\n")
