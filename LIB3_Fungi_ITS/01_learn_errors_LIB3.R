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

cat("Number of samples:", length(sample.names), "\n")

dir.create("dada2_out", showWarnings = FALSE)

# Learn forward error rates
errF <- learnErrors(
  fnFs,
  multithread = TRUE
)

# Learn reverse error rates
errR <- learnErrors(
  fnRs,
  multithread = TRUE
)

saveRDS(errF, "dada2_out/error_model_F.rds")
saveRDS(errR, "dada2_out/error_model_R.rds")

pdf("dada2_out/error_rates_F.pdf")
plotErrors(errF, nominalQ = TRUE)
dev.off()

pdf("dada2_out/error_rates_R.pdf")
plotErrors(errR, nominalQ = TRUE)
dev.off()

cat("\nError-rate learning completed.\n")
