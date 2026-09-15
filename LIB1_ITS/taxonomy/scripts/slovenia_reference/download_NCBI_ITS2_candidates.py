#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path
import csv
import time

# --------------------------------------------------
# SETTINGS
# --------------------------------------------------

Entrez.email = "your_email@example.com"

# Ako nemaš NCBI API key, ostavi None
Entrez.api_key = None

input_file = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "Slovenia_missing_taxa_NCBI_ITS2_hits.tsv"
)

output_fasta = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "NCBI_ITS2_candidates_raw.fasta"
)

output_metadata = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "NCBI_ITS2_candidates_metadata.tsv"
)

log_file = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "NCBI_ITS2_download.log"
)

# Conservative delay
delay = 0.4 if Entrez.api_key is None else 0.12


# --------------------------------------------------
# READ TAXA WITH HITS
# --------------------------------------------------

species_list = []

with open(input_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        species = row["Species"].strip()
        hits = int(row["NCBI_ITS2_hits"])

        if hits > 0:
            species_list.append(species)

print(f"Taxa with NCBI ITS2 hits: {len(species_list)}")


# --------------------------------------------------
# DOWNLOAD
# --------------------------------------------------

metadata_rows = []
downloaded = 0
failed_species = 0

with open(output_fasta, "w") as fasta_out, \
     open(log_file, "w") as log:

    for i, species in enumerate(species_list, start=1):

        query = (
            f'"{species}"[Organism] AND '
            f'("internal transcribed spacer 2"[Title] OR ITS2[Title])'
        )

        print(f"[{i}/{len(species_list)}] {species}", flush=True)

        try:

            # Search accession IDs
            handle = Entrez.esearch(
                db="nuccore",
                term=query,
                retmax=10000
            )

            search_result = Entrez.read(handle)
            handle.close()

            ids = search_result["IdList"]

            if not ids:
                continue

            time.sleep(delay)

            # Fetch GenBank records so that we retain metadata
            # Do in batches to avoid huge requests
            batch_size = 100

            for start in range(0, len(ids), batch_size):

                batch_ids = ids[start:start + batch_size]

                handle = Entrez.efetch(
                    db="nuccore",
                    id=",".join(batch_ids),
                    rettype="gb",
                    retmode="text"
                )

                records = list(
                    SeqIO.parse(handle, "genbank")
                )

                handle.close()

                for record in records:

                    accession = record.id
                    description = record.description
                    length = len(record.seq)

                    organism = ""
                    taxonomy = ""

                    if "organism" in record.annotations:
                        organism = record.annotations["organism"]

                    if "taxonomy" in record.annotations:
                        taxonomy = ";".join(
                            record.annotations["taxonomy"]
                        )

                    # Preserve original sequence exactly
                    SeqIO.write(
                        record,
                        fasta_out,
                        "fasta"
                    )

                    metadata_rows.append({
                        "Requested_species": species,
                        "Accession": accession,
                        "Organism": organism,
                        "Length": length,
                        "Description": description,
                        "NCBI_taxonomy": taxonomy
                    })

                    downloaded += 1

                time.sleep(delay)

            log.write(
                f"{species}\t{len(ids)}\tOK\n"
            )
            log.flush()

        except Exception as e:

            failed_species += 1

            log.write(
                f"{species}\tERROR\t{e}\n"
            )
            log.flush()

            print(f"  ERROR: {e}")

            time.sleep(2)


# --------------------------------------------------
# SAVE METADATA
# --------------------------------------------------

with open(output_metadata, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Requested_species",
            "Accession",
            "Organism",
            "Length",
            "Description",
            "NCBI_taxonomy"
        ],
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(metadata_rows)


print()
print("=== NCBI ITS2 DOWNLOAD COMPLETE ===")
print(f"Taxa attempted:        {len(species_list)}")
print(f"Records downloaded:    {downloaded}")
print(f"Species with errors:   {failed_species}")
print()
print(f"Raw FASTA: {output_fasta}")
print(f"Metadata:  {output_metadata}")
print(f"Log:       {log_file}")
