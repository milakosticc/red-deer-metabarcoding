library(dada2)

# Input folder: quality-filtered fungal reads
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

# Extract sample names
sample.names <- sub("_filtered_R1.fastq.gz$", "", basename(fnFs))

# Check that R1/R2 files correspond
stopifnot(length(fnFs) == length(fnRs))

sample.names.R <- sub("_filtered_R2.fastq.gz$", "", basename(fnRs))
stopifnot(all(sample.names == sample.names.R))

# TEST ONLY: use first 4 samples
fnFs.test <- fnFs[1:4]
fnRs.test <- fnRs[1:4]
samples.test <- sample.names[1:4]

names(fnFs.test) <- samples.test
names(fnRs.test) <- samples.test

cat("Testing samples:\n")
print(samples.test)

# Learn error rates
errF <- learnErrors(
  fnFs.test,
  multithread = TRUE
)

errR <- learnErrors(
  fnRs.test,
  multithread = TRUE
)

saveRDS(errF, "test_dada2_out/error_model_F.rds")
saveRDS(errR, "test_dada2_out/error_model_R.rds")

# Dereplication
derepFs <- derepFastq(fnFs.test)
derepRs <- derepFastq(fnRs.test)

names(derepFs) <- samples.test
names(derepRs) <- samples.test

# Denoising
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

# Merge paired reads
mergers <- mergePairs(
  dadaFs,
  derepFs,
  dadaRs,
  derepRs,
  verbose = TRUE
)

# ASV sequence table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# Tracking table
getN <- function(x) sum(getUniques(x))

track <- cbind(
  denoisedF = sapply(dadaFs, getN),
  denoisedR = sapply(dadaRs, getN),
  merged = sapply(mergers, getN),
  nonchim = rowSums(seqtab.nochim)
)

write.csv(
  track,
  "test_dada2_out/test_dada2_tracking.csv",
  quote = FALSE
)

# Save test sequence table
saveRDS(
  seqtab.nochim,
  "test_dada2_out/test_seqtab_nochim.rds"
)

# FASTA for ITSx test
asv.seqs <- colnames(seqtab.nochim)

fasta <- character(length(asv.seqs) * 2)

fasta[seq(1, length(fasta), 2)] <-
  paste0(">ASV", seq_along(asv.seqs))

fasta[seq(2, length(fasta), 2)] <- asv.seqs

writeLines(
  fasta,
  "test_dada2_out/test_ASVs.fasta"
)

cat("\nTEST COMPLETED\n")
cat("Sequence table dimensions:\n")
print(dim(seqtab))

cat("Non-chimeric sequence table dimensions:\n")
print(dim(seqtab.nochim))

cat("\nTracking table:\n")
print(track)
