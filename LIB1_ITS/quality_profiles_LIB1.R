library(dada2)

path <- "working/Plants_ITS/primersCut_out"

fnFs <- sort(list.files(
  path,
  pattern = "_trimmed_R1.fastq.gz$",
  full.names = TRUE
))

fnRs <- sort(list.files(
  path,
  pattern = "_trimmed_R2.fastq.gz$",
  full.names = TRUE
))

cat("Number of R1 files:", length(fnFs), "\n")
cat("Number of R2 files:", length(fnRs), "\n")

if (length(fnFs) != length(fnRs)) {
  stop("Different number of R1 and R2 files!")
}

pdf(
  "working/Plants_ITS/quality_profiles_LIB1.pdf",
  width = 12,
  height = 8
)

print(plotQualityProfile(fnFs))
print(plotQualityProfile(fnRs))

dev.off()
