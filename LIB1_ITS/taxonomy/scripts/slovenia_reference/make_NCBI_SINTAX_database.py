#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from pathlib import Path
from collections import Counter, defaultdict
import csv
import re
import time


# ============================================================
# SETTINGS
# ============================================================

Entrez.email = "your_email@example.com"

base = Path("taxonomy/reference/Slovenia_ITS2")
ncbi = base / "NCBI"

bold_fasta = base / "Slovenia_BOLD_Plantae_ITS2_SINTAX.fasta"

validated_fasta = ncbi / "NCBI_ITS2_validated.fasta"
validated_metadata = ncbi / "NCBI_ITS2_validated.tsv"

output_fasta = ncbi / "NCBI_Plantae_ITS2_SINTAX.fasta"
output_mapping = ncbi / "NCBI_Plantae_ITS2_SINTAX_mapping.tsv"
overlap_file = ncbi / "NCBI_excluded_already_in_BOLD.tsv"

delay = 0.4


# ============================================================
# HELPERS
# ============================================================

def clean_taxon(x):

    if not x:
        return ""

    x = str(x).strip()

    x = re.sub(r"\s+", "_", x)
    x = x.replace(";", "")
    x = x.replace(",", "")

    return x


def parse_bold_taxonomy(header):

    if "tax=" not in header:
        return {}

    tax = header.split("tax=", 1)[1]

    result = {}

    for item in tax.split(","):

        item = item.rstrip(";")

        if ":" not in item:
            continue

        rank, value = item.split(":", 1)

        result[rank] = value

    return result


# ============================================================
# READ BOLD TAXONOMY
# ============================================================

bold_species = set()

genus_lineages = defaultdict(list)

for record in SeqIO.parse(bold_fasta, "fasta"):

    tax = parse_bold_taxonomy(record.description)

    species = tax.get("s", "").replace("_", " ")
    genus = tax.get("g", "")

    if species:
        bold_species.add(species.lower())

    if genus:

        lineage = (
            tax.get("p", ""),
            tax.get("c", ""),
            tax.get("o", ""),
            tax.get("f", "")
        )

        genus_lineages[genus].append(lineage)


# Consensus BOLD lineage for each genus

bold_genus_consensus = {}

for genus, lineages in genus_lineages.items():

    bold_genus_consensus[genus] = (
        Counter(lineages).most_common(1)[0][0]
    )


print(f"BOLD species represented: {len(bold_species)}")
print(f"BOLD genera represented:  {len(bold_genus_consensus)}")


# ============================================================
# READ VALIDATED NCBI METADATA
# ============================================================

metadata = {}

taxids = set()

with open(validated_metadata, newline="") as f:

    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:

        accession = row["Accession"].strip()

        metadata[accession] = row

        taxid = row["Record_TaxID"].strip()

        if taxid:
            taxids.add(taxid)


print(f"Validated NCBI records: {len(metadata)}")
print(f"Unique NCBI TaxIDs:     {len(taxids)}")


# ============================================================
# READ SEQUENCES
# ============================================================

sequences = {}

for record in SeqIO.parse(validated_fasta, "fasta"):

    accession = record.id.split("|")[0]

    sequences[accession] = str(record.seq).upper()


# ============================================================
# FETCH NCBI TAXONOMY
# ============================================================

taxonomy_by_taxid = {}

taxid_list = sorted(taxids)

batch_size = 100


for start in range(0, len(taxid_list), batch_size):

    batch = taxid_list[start:start + batch_size]

    print(
        f"Fetching taxonomy "
        f"{start + 1}-"
        f"{min(start + batch_size, len(taxid_list))}"
        f"/{len(taxid_list)}",
        flush=True
    )

    handle = Entrez.efetch(
        db="taxonomy",
        id=",".join(batch),
        retmode="xml"
    )

    records = Entrez.read(handle)
    handle.close()

    for taxon in records:

        taxid = str(taxon["TaxId"])

        rank_map = {}

        for node in taxon.get("LineageEx", []):

            rank = str(node.get("Rank", ""))
            name = str(node.get("ScientificName", ""))

            if rank and name:
                rank_map[rank] = name

        current_rank = str(
            taxon.get("Rank", "")
        )

        current_name = str(
            taxon.get("ScientificName", "")
        )

        if current_rank and current_name:
            rank_map[current_rank] = current_name

        taxonomy_by_taxid[taxid] = rank_map

    time.sleep(delay)


# ============================================================
# CONVERT TO SINTAX
# ============================================================

accepted = []
excluded_overlap = []


for accession, row in metadata.items():

    taxid = row["Record_TaxID"].strip()

    rank_map = taxonomy_by_taxid.get(
        taxid,
        {}
    )

    # Accepted NCBI species
    species = rank_map.get("species", "")

    # Conservative fallback
    if not species:

        organism = row["Organism"].strip()
        parts = organism.split()

        if len(parts) >= 2:
            species = " ".join(parts[:2])

    genus = rank_map.get("genus", "")

    if not genus and species:
        genus = species.split()[0]


    # --------------------------------------------------------
    # Do not add NCBI taxon if accepted species
    # already exists in BOLD
    # --------------------------------------------------------

    if species.lower() in bold_species:

        excluded_overlap.append({
            "Accession": accession,
            "Requested_species": row["Requested_species"],
            "NCBI_accepted_species": species,
            "Reason": "Accepted species already represented in BOLD"
        })

        continue


    # --------------------------------------------------------
    # Prefer BOLD higher taxonomy for the genus
    # so the combined database uses the same backbone
    # --------------------------------------------------------

    genus_key = genus.replace(" ", "_")

    if genus_key in bold_genus_consensus:

        phylum, class_name, order, family = (
            bold_genus_consensus[genus_key]
        )

        taxonomy_source = "BOLD_GENUS_BACKBONE"

    else:

        phylum = rank_map.get("phylum", "")
        class_name = rank_map.get("class", "")
        order = rank_map.get("order", "")
        family = rank_map.get("family", "")

        taxonomy_source = "NCBI_LINEAGE"


    if not species or not genus:

        continue


    sequence = sequences.get(accession)

    if not sequence:
        continue


    accepted.append({
        "Accession": accession,
        "Requested_species": row["Requested_species"],
        "NCBI_organism": row["Organism"],
        "NCBI_accepted_species": species,
        "Phylum": phylum,
        "Class": class_name,
        "Order": order,
        "Family": family,
        "Genus": genus,
        "Taxonomy_source": taxonomy_source,
        "Sequence": sequence
    })


# ============================================================
# WRITE SINTAX FASTA
# ============================================================

with open(output_fasta, "w") as f:

    for i, row in enumerate(
        accepted,
        start=1
    ):

        ref_id = f"NCBI_ITS2_REF{i:06d}"

        header = (
            f">{ref_id};tax="
            f"k:Plantae,"
            f"p:{clean_taxon(row['Phylum'])},"
            f"c:{clean_taxon(row['Class'])},"
            f"o:{clean_taxon(row['Order'])},"
            f"f:{clean_taxon(row['Family'])},"
            f"g:{clean_taxon(row['Genus'])},"
            f"s:{clean_taxon(row['NCBI_accepted_species'])};"
        )

        f.write(header + "\n")

        seq = row["Sequence"]

        for j in range(0, len(seq), 80):
            f.write(seq[j:j+80] + "\n")


# ============================================================
# WRITE MAPPING
# ============================================================

with open(output_mapping, "w", newline="") as f:

    fields = [
        "Accession",
        "Requested_species",
        "NCBI_organism",
        "NCBI_accepted_species",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Taxonomy_source"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()

    for row in accepted:

        writer.writerow({
            key: row[key]
            for key in fields
        })


# ============================================================
# WRITE BOLD OVERLAPS
# ============================================================

with open(overlap_file, "w", newline="") as f:

    fields = [
        "Accession",
        "Requested_species",
        "NCBI_accepted_species",
        "Reason"
    ]

    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(excluded_overlap)


# ============================================================
# SUMMARY
# ============================================================

accepted_species = {
    row["NCBI_accepted_species"]
    for row in accepted
}

print()
print("=== NCBI SINTAX DATABASE COMPLETE ===")
print(f"Validated NCBI sequences:       {len(metadata)}")
print(f"Excluded because already BOLD: {len(excluded_overlap)}")
print(f"NCBI sequences added:           {len(accepted)}")
print(f"Additional accepted species:    {len(accepted_species)}")
print()
print(f"FASTA:   {output_fasta}")
print(f"Mapping: {output_mapping}")
print(f"Overlap: {overlap_file}")
