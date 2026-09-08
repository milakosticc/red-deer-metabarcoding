library(dada2)

input_dir <- "working/Fungi_ITS/primersCut_out"
output_dir <- "working/Fungi_ITS/qualFiltered_out"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

fnFs <- sort(list.files(
  input_dir,
  pattern = "_trimmed_R1.fastq.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  input_dir,
  pattern = "_trimmed_R2.fastq.gz$",
  full.names = TRUE
))

sampleF <- sub(
  "_trimmed_R1.fastq.gz$",
  "",
  basename(fnFs)
)

sampleR <- sub(
  "_trimmed_R2.fastq.gz$",
  "",
  basename(fnRs)
)

if (!identical(sampleF, sampleR)) {
  stop("R1 and R2 sample names do not match!")
}

cat("Number of samples:", length(sampleF), "\n")

filtFs <- file.path(
  output_dir,
  paste0(sampleF, "_filtered_R1.fastq.gz")
)

filtRs <- file.path(
  output_dir,
  paste0(sampleF, "_filtered_R2.fastq.gz")
)

out <- filterAndTrim(
  fnFs,
  filtFs,
  fnRs,
  filtRs,
  maxN = 0,
  maxEE = c(2,2),
  truncQ = 2,
  truncLen = c(0,0),
  minLen = 100,
  maxLen = 600,
  matchIDs = TRUE,
  compress = TRUE,
  multithread = TRUE
)

summary <- data.frame(
  sample = sampleF,
  reads_in = out[, "reads.in"],
  reads_out = out[, "reads.out"]
)

summary$retained_percent <- round(
  100 * summary$reads_out / summary$reads_in,
  2
)

write.csv(
  summary,
  "working/Fungi_ITS/seq_count_summary.csv",
  row.names = FALSE
)

print(summary)
