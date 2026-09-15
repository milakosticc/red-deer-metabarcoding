#!/usr/bin/env python3

from Bio import SeqIO
from pathlib import Path
import csv


base = Path("taxonomy/reference/Slovenia_ITS2/NCBI")

input_fasta = base / "NCBI_ITS2_extracted_v2_preQC.fasta"
input_metadata = base / "NCBI_ITS2_extracted_v2_preQC.tsv"

output_fasta = base / "NCBI_ITS2_validated.fasta"
output_metadata = base / "NCBI_ITS2_validated.tsv"

MIN_LEN = 100
MAX_LEN = 700


# --------------------------------------------------
# Read metadata
# --------------------------------------------------

metadata = {}

with open(input_metadata, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        metadata[row["Accession"]] = row


# --------------------------------------------------
# Filter FASTA
# --------------------------------------------------

kept = []
rejected_length = 0

for record in SeqIO.parse(input_fasta, "fasta"):

    accession = record.id.split("|")[0]

    length = len(record.seq)

    if length < MIN_LEN or length > MAX_LEN:
        rejected_length += 1
        continue

    kept.append(record)


SeqIO.write(
    kept,
    output_fasta,
    "fasta"
)


# --------------------------------------------------
# Save matching metadata
# --------------------------------------------------

kept_accessions = {
    record.id.split("|")[0]
    for record in kept
}

fields = [
    "Requested_species",
    "Accession",
    "Organism",
    "Record_TaxID",
    "Taxonomy_match",
    "Extraction_method",
    "Feature_type",
    "Length",
    "Description"
]

kept_rows = []

for accession in kept_accessions:

    if accession in metadata:
        kept_rows.append(metadata[accession])


kept_rows.sort(
    key=lambda x: (
        x["Requested_species"],
        x["Accession"]
    )
)


with open(output_metadata, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(kept_rows)


taxa = {
    row["Requested_species"]
    for row in kept_rows
}


print("=== NCBI ITS2 VALIDATION COMPLETE ===")
print(f"Input sequences:       {len(metadata)}")
print(f"Validated sequences:   {len(kept)}")
print(f"Removed by length:     {rejected_length}")
print(f"Taxa represented:      {len(taxa)}")
print()
print(f"FASTA:    {output_fasta}")
print(f"Metadata: {output_metadata}")
