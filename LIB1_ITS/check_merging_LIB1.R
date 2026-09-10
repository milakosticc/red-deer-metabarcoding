library(dada2)

base <- "/home/mila/aapraksa/samples/working/Plants_ITS"

filtered <- file.path(base, "qualFiltered_out")
denoised <- file.path(base, "dada2_out/denoised")

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

dadaFs <- readRDS(file.path(denoised, "LIB1_dada_R1.rds"))
dadaRs <- readRDS(file.path(denoised, "LIB1_dada_R2.rds"))

mergers.all <- mergePairs(
  dadaFs,
  fnFs,
  dadaRs,
  fnRs,
  returnRejects = TRUE
)

names(mergers.all) <- sample.names

for (s in c("LME1676", "LME1685", "LME1690", "LME1691")) {
  cat("\n===== ", s, " =====\n", sep="")
  print(
    head(
      mergers.all[[s]][
        order(-mergers.all[[s]]$abundance),
      ],
      15
    )
  )
}
