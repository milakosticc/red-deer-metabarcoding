#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path
import csv
import re
import time


# ============================================================
# SETTINGS
# ============================================================

Entrez.email = "your_email@example.com"
Entrez.api_key = None

delay = 0.38 if Entrez.api_key is None else 0.12


# ============================================================
# INPUT
# ============================================================

metadata_file = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "NCBI_ITS2_candidates_metadata.tsv"
)


# ============================================================
# OUTPUT
# ============================================================

outdir = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI"
)

output_fasta = outdir / "NCBI_ITS2_extracted_v2_preQC.fasta"

output_metadata = outdir / "NCBI_ITS2_extracted_v2_preQC.tsv"

rejected_file = outdir / "NCBI_ITS2_v2_rejected.tsv"

taxonomy_cache_file = outdir / "NCBI_requested_species_taxids.tsv"


# ============================================================
# HELPERS
# ============================================================

def normalize_binomial(name):

    parts = name.strip().split()

    if len(parts) >= 2:
        return " ".join(parts[:2]).lower()

    return name.strip().lower()


def get_record_taxid(record):

    for feature in record.features:

        if feature.type != "source":
            continue

        for value in feature.qualifiers.get("db_xref", []):

            if value.startswith("taxon:"):
                return value.split(":", 1)[1]

    return ""


def is_exact_its2_feature(feature):
    """
    Accept only a feature whose PRODUCT explicitly describes ITS2.

    We intentionally DO NOT use 'note' because a broad feature can
    say that it contains ITS1 + 5.8S + ITS2 + other regions.
    """

    products = feature.qualifiers.get("product", [])

    if not products:
        return False

    text = " ".join(str(x) for x in products).lower()

    has_its2 = (
        "internal transcribed spacer 2" in text
        or re.search(r"\bits[\s_-]*2\b", text)
    )

    has_its1 = (
        "internal transcribed spacer 1" in text
        or re.search(r"\bits[\s_-]*1\b", text)
    )

    has_58s = "5.8s" in text

    return bool(
        has_its2
        and not has_its1
        and not has_58s
    )


def whole_record_is_its2_only(description):

    text = description.lower()

    has_its2 = (
        "internal transcribed spacer 2" in text
        or re.search(r"\bits[\s_-]*2\b", text)
    )

    has_its1 = (
        "internal transcribed spacer 1" in text
        or re.search(r"\bits[\s_-]*1\b", text)
    )

    other_regions = any([
        "5.8s" in text,
        "small subunit ribosomal" in text,
        "large subunit ribosomal" in text,
        "internal transcribed spacer 1" in text,
    ])

    return bool(
        has_its2
        and not has_its1
        and not other_regions
    )


# ============================================================
# READ METADATA
# ============================================================

metadata = {}

requested_species_set = set()

with open(metadata_file, newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        accession = row["Accession"].strip()

        if not accession:
            continue

        # Duplicate accession IDs are intentionally collapsed.
        metadata[accession] = row

        requested_species_set.add(
            row["Requested_species"].strip()
        )


accessions = list(metadata.keys())

requested_species_list = sorted(requested_species_set)

print(f"Unique accessions: {len(accessions)}")
print(f"Requested taxa:    {len(requested_species_list)}")


# ============================================================
# BUILD / READ NCBI TAXONOMY SYNONYM MAP
#
# We ask NCBI taxonomy which TaxIDs correspond to the name in
# the Slovenian whitelist. This allows:
#
# Acer ginnala
#       ->
# Acer tataricum subsp. ginnala
#
# when NCBI recognizes Acer ginnala as a synonym.
# ============================================================

species_taxids = {}


if taxonomy_cache_file.exists():

    print()
    print("Reading existing taxonomy cache...")

    with open(taxonomy_cache_file, newline="") as f:

        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:

            species = row["Requested_species"]

            taxids = {
                x
                for x in row["NCBI_TaxIDs"].split(",")
                if x
            }

            species_taxids[species] = taxids


else:

    print()
    print("Building NCBI synonym/taxonomy map...")

    cache_rows = []

    for i, species in enumerate(
        requested_species_list,
        start=1
    ):

        print(
            f"[taxonomy {i}/{len(requested_species_list)}] "
            f"{species}",
            flush=True
        )

        taxids = set()

        try:

            handle = Entrez.esearch(
                db="taxonomy",
                term=f'"{species}"[All Names]',
                retmax=20
            )

            result = Entrez.read(handle)
            handle.close()

            taxids = set(result["IdList"])

        except Exception as e:

            print(f"  TAXONOMY ERROR: {e}")

        species_taxids[species] = taxids

        cache_rows.append({
            "Requested_species": species,
            "NCBI_TaxIDs": ",".join(sorted(taxids))
        })

        time.sleep(delay)


    with open(
        taxonomy_cache_file,
        "w",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            fieldnames=[
                "Requested_species",
                "NCBI_TaxIDs"
            ],
            delimiter="\t",
            lineterminator="\n"
        )

        writer.writeheader()
        writer.writerows(cache_rows)


# ============================================================
# PROCESS NCBI RECORDS
# ============================================================

accepted = []
rejected = []

batch_size = 100


for start in range(
    0,
    len(accessions),
    batch_size
):

    batch = accessions[
        start:start + batch_size
    ]

    first = start + 1
    last = min(
        start + batch_size,
        len(accessions)
    )

    print(
        f"Processing {first}-{last}/{len(accessions)}",
        flush=True
    )


    try:

        handle = Entrez.efetch(
            db="nuccore",
            id=",".join(batch),
            rettype="gb",
            retmode="text"
        )

        records = list(
            SeqIO.parse(handle, "genbank")
        )

        handle.close()


    except Exception as e:

        for accession in batch:

            rejected.append({
                "Requested_species":
                    metadata[accession]["Requested_species"],
                "Accession": accession,
                "Organism": "",
                "Record_TaxID": "",
                "Reason": f"FETCH_ERROR: {e}"
            })

        time.sleep(2)
        continue


    returned = set()


    for record in records:

        accession = record.id
        returned.add(accession)

        info = metadata.get(accession)

        if info is None:
            continue


        requested_species = (
            info["Requested_species"].strip()
        )

        organism = (
            record.annotations
            .get("organism", "")
            .strip()
        )

        record_taxid = get_record_taxid(
            record
        )


        # ====================================================
        # TAXONOMIC VALIDATION
        # ====================================================

        requested_taxids = species_taxids.get(
            requested_species,
            set()
        )

        exact_binomial_match = (
            normalize_binomial(requested_species)
            ==
            normalize_binomial(organism)
        )

        synonym_taxid_match = (
            bool(record_taxid)
            and record_taxid in requested_taxids
        )


        if not (
            exact_binomial_match
            or synonym_taxid_match
        ):

            rejected.append({
                "Requested_species": requested_species,
                "Accession": accession,
                "Organism": organism,
                "Record_TaxID": record_taxid,
                "Reason": "TAXONOMY_MISMATCH"
            })

            continue


        # ====================================================
        # EXACT ITS2 FEATURE
        # ====================================================

        its2_features = []

        for feature in record.features:

            if feature.type == "source":
                continue

            if not is_exact_its2_feature(feature):
                continue

            try:

                seq = feature.extract(
                    record.seq
                )

            except Exception:
                continue

            if len(seq) == 0:
                continue

            its2_features.append(
                (
                    len(seq),
                    feature.type,
                    seq
                )
            )


        extracted_seq = None
        extraction_method = ""
        feature_type = ""


        if its2_features:

            # Normally only one exact ITS2 feature exists.
            # If several exist, use the shortest exact feature.
            its2_features.sort(
                key=lambda x: x[0]
            )

            (
                feature_length,
                feature_type,
                extracted_seq
            ) = its2_features[0]

            extraction_method = (
                "EXACT_ANNOTATED_ITS2"
            )


        # ====================================================
        # ITS2-ONLY WHOLE RECORD FALLBACK
        # ====================================================

        elif whole_record_is_its2_only(
            record.description
        ):

            if len(record.seq) > 0:

                extracted_seq = record.seq

                feature_type = "whole_record"

                extraction_method = (
                    "ITS2_ONLY_RECORD"
                )


        # ====================================================
        # NO SAFE ITS2 REGION
        # ====================================================

        if extracted_seq is None:

            rejected.append({
                "Requested_species": requested_species,
                "Accession": accession,
                "Organism": organism,
                "Record_TaxID": record_taxid,
                "Reason": "NO_EXACT_ITS2_REGION"
            })

            continue


        sequence = str(
            extracted_seq
        ).upper()


        if not sequence:

            rejected.append({
                "Requested_species": requested_species,
                "Accession": accession,
                "Organism": organism,
                "Record_TaxID": record_taxid,
                "Reason": "EMPTY_SEQUENCE"
            })

            continue


        accepted.append({
            "Requested_species": requested_species,
            "Accession": accession,
            "Organism": organism,
            "Record_TaxID": record_taxid,
            "Taxonomy_match":
                (
                    "NCBI_SYNONYM"
                    if synonym_taxid_match
                    and not exact_binomial_match
                    else "NAME_MATCH"
                ),
            "Extraction_method":
                extraction_method,
            "Feature_type":
                feature_type,
            "Length":
                len(sequence),
            "Description":
                record.description,
            "Sequence":
                sequence
        })


    # ========================================================
    # ACCESSIONS NOT RETURNED
    # ========================================================

    for accession in batch:

        if accession not in returned:

            rejected.append({
                "Requested_species":
                    metadata[accession]["Requested_species"],
                "Accession":
                    accession,
                "Organism":
                    "",
                "Record_TaxID":
                    "",
                "Reason":
                    "NOT_RETURNED_BY_NCBI"
            })


    time.sleep(delay)


# ============================================================
# WRITE FASTA
# ============================================================

with open(output_fasta, "w") as f:

    for row in accepted:

        species_header = (
            row["Requested_species"]
            .replace(" ", "_")
        )

        f.write(
            f">{row['Accession']}"
            f"|species={species_header}"
            f"|source=NCBI"
            f"|method={row['Extraction_method']}\n"
        )

        seq = row["Sequence"]

        for i in range(
            0,
            len(seq),
            80
        ):

            f.write(
                seq[i:i + 80] + "\n"
            )


# ============================================================
# WRITE METADATA
# ============================================================

metadata_fields = [
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


with open(
    output_metadata,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=metadata_fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in accepted:

        writer.writerow({
            field: row[field]
            for field in metadata_fields
        })


# ============================================================
# WRITE REJECTED
# ============================================================

with open(
    rejected_file,
    "w",
    newline=""
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Requested_species",
            "Accession",
            "Organism",
            "Record_TaxID",
            "Reason"
        ],
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(rejected)


# ============================================================
# SUMMARY
# ============================================================

represented_taxa = {
    row["Requested_species"]
    for row in accepted
}

synonym_records = sum(
    1 for row in accepted
    if row["Taxonomy_match"] == "NCBI_SYNONYM"
)

name_match_records = sum(
    1 for row in accepted
    if row["Taxonomy_match"] == "NAME_MATCH"
)


print()
print("=== NCBI ITS2 V2 EXTRACTION COMPLETE ===")

print(f"Unique accessions processed: {len(accessions)}")

print(
    f"ITS2 sequences extracted:   {len(accepted)}"
)

print(
    f"Taxa represented:           {len(represented_taxa)}"
)

print(
    f"Direct name matches:        {name_match_records}"
)

print(
    f"NCBI synonym matches:       {synonym_records}"
)

print(
    f"Records rejected:           {len(rejected)}"
)

print()
print(f"FASTA:    {output_fasta}")
print(f"Metadata: {output_metadata}")
print(f"Rejected: {rejected_file}")
