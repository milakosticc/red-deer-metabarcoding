#!/usr/bin/env python3

import csv
from collections import defaultdict
from pathlib import Path

taxonomy_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_OTU97_Slovenia_taxonomy.tsv"
)

abundance_file = Path(
    "taxonomy/otu_97/LIB1_OTU97_abundance_counts.tsv"
)

top5_output = Path(
    "taxonomy/otu_97/slovenia/LIB1_Slovenia_top5_taxa_by_sample.tsv"
)


# -----------------------------
# Read taxonomy
# -----------------------------

taxonomy = {}

with open(taxonomy_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        taxonomy[row["OTU"]] = row

print(f"Taxonomy OTUs: {len(taxonomy)}")


# -----------------------------
# Read abundance table
# FORMAT: samples in rows, OTUs in columns
# -----------------------------

with open(abundance_file, newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    columns = reader.fieldnames

sample_column = columns[0]
otu_columns = columns[1:]

print(f"Samples: {len(rows)}")
print(f"OTU columns: {len(otu_columns)}")
print(f"Sample column: {sample_column}")


# -----------------------------
# Check OTU matching
# -----------------------------

otu_set = set(otu_columns)
taxonomy_set = set(taxonomy)

print(f"OTUs without taxonomy: {len(otu_set - taxonomy_set)}")
print(f"Taxonomy OTUs absent from abundance table: {len(taxonomy_set - otu_set)}")


# -----------------------------
# Helpers
# -----------------------------

def clean_taxon(value):

    if value is None:
        return ""

    value = value.strip()

    if value.lower() in {
        "",
        "na",
        "nan",
        "none",
        "unclassified",
        "unassigned"
    }:
        return ""

    return value


def format_result(taxon, count, total):

    percent = (count / total * 100) if total > 0 else 0

    return f"{taxon} ({count:,}; {percent:.2f}%)"


# -----------------------------
# Build Top 5 table
# -----------------------------

summary_rows = []

for row in rows:

    sample = row[sample_column]

    species_counts = defaultdict(int)
    genus_counts = defaultdict(int)
    family_counts = defaultdict(int)

    total_reads = 0

    for otu in otu_columns:

        try:
            count = int(float(row[otu]))
        except (ValueError, TypeError):
            count = 0

        total_reads += count

        if count == 0:
            continue

        tax = taxonomy.get(otu)

        if not tax:
            continue

        species = clean_taxon(tax.get("Species", ""))
        genus = clean_taxon(tax.get("Genus", ""))
        family = clean_taxon(tax.get("Family", ""))

        if species:
            species_counts[species] += count

        if genus:
            genus_counts[genus] += count

        if family:
            family_counts[family] += count


    def get_top5(counts):

        ordered = sorted(
            counts.items(),
            key=lambda x: (-x[1], x[0])
        )[:5]

        results = [
            format_result(taxon, count, total_reads)
            for taxon, count in ordered
        ]

        while len(results) < 5:
            results.append("")

        return results


    top_species = get_top5(species_counts)
    top_genus = get_top5(genus_counts)
    top_family = get_top5(family_counts)

    result = {
        "Sample": sample
    }

    for i in range(5):
        result[f"Top_{i+1}_Species"] = top_species[i]

    for i in range(5):
        result[f"Top_{i+1}_Genus"] = top_genus[i]

    for i in range(5):
        result[f"Top_{i+1}_Family"] = top_family[i]

    summary_rows.append(result)


# -----------------------------
# Save
# -----------------------------

columns = (
    ["Sample"]
    + [f"Top_{i}_Species" for i in range(1, 6)]
    + [f"Top_{i}_Genus" for i in range(1, 6)]
    + [f"Top_{i}_Family" for i in range(1, 6)]
)

with open(top5_output, "w", newline="") as out:

    writer = csv.DictWriter(
        out,
        fieldnames=columns,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(summary_rows)


print()
print("Created:")
print(top5_output)
