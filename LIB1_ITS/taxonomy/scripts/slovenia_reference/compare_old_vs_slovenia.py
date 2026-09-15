#!/usr/bin/env python3

import csv
import re
from pathlib import Path

old_file = Path(
    "taxonomy/final/LIB1_top5_species_genus_family_per_sample.tsv"
)

new_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_Slovenia_top5_taxa_by_sample.tsv"
)

output_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_old_vs_Slovenia_comparison.tsv"
)


def read_table(path):
    data = {}

    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")

        for row in reader:
            data[row["Sample"]] = row

    return data


def parse_taxon(value):

    if not value:
        return "", None

    # Example:
    # Rubus (22,980; 80.03%)

    match = re.match(
        r"^(.*?)\s+\([\d,]+;\s*([\d.]+)%\)$",
        value.strip()
    )

    if not match:
        return value.strip(), None

    taxon = match.group(1).strip()
    percent = float(match.group(2))

    return taxon, percent


old = read_table(old_file)
new = read_table(new_file)

samples = sorted(set(old) & set(new))

rows = []

for sample in samples:

    old_species, old_species_pct = parse_taxon(
        old[sample].get("Top_1_Species", "")
    )

    new_species, new_species_pct = parse_taxon(
        new[sample].get("Top_1_Species", "")
    )

    old_genus, old_genus_pct = parse_taxon(
        old[sample].get("Top_1_Genus", "")
    )

    new_genus, new_genus_pct = parse_taxon(
        new[sample].get("Top_1_Genus", "")
    )

    old_family, old_family_pct = parse_taxon(
        old[sample].get("Top_1_Family", "")
    )

    new_family, new_family_pct = parse_taxon(
        new[sample].get("Top_1_Family", "")
    )


    genus_diff = ""
    if old_genus_pct is not None and new_genus_pct is not None:
        genus_diff = round(new_genus_pct - old_genus_pct, 2)

    family_diff = ""
    if old_family_pct is not None and new_family_pct is not None:
        family_diff = round(new_family_pct - old_family_pct, 2)


    rows.append({
        "Sample": sample,

        "Old_Top_Species": old_species,
        "Slovenia_Top_Species": new_species,
        "Same_Top_Species": old_species == new_species,

        "Old_Top_Genus": old_genus,
        "Old_Genus_percent": old_genus_pct,
        "Slovenia_Top_Genus": new_genus,
        "Slovenia_Genus_percent": new_genus_pct,
        "Same_Top_Genus": old_genus == new_genus,
        "Genus_percent_difference": genus_diff,

        "Old_Top_Family": old_family,
        "Old_Family_percent": old_family_pct,
        "Slovenia_Top_Family": new_family,
        "Slovenia_Family_percent": new_family_pct,
        "Same_Top_Family": old_family == new_family,
        "Family_percent_difference": family_diff
    })


fieldnames = [
    "Sample",

    "Old_Top_Species",
    "Slovenia_Top_Species",
    "Same_Top_Species",

    "Old_Top_Genus",
    "Old_Genus_percent",
    "Slovenia_Top_Genus",
    "Slovenia_Genus_percent",
    "Same_Top_Genus",
    "Genus_percent_difference",

    "Old_Top_Family",
    "Old_Family_percent",
    "Slovenia_Top_Family",
    "Slovenia_Family_percent",
    "Same_Top_Family",
    "Family_percent_difference"
]


with open(output_file, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=fieldnames,
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(rows)


same_species = sum(
    1 for r in rows
    if r["Same_Top_Species"]
)

same_genus = sum(
    1 for r in rows
    if r["Same_Top_Genus"]
)

same_family = sum(
    1 for r in rows
    if r["Same_Top_Family"]
)


print("=== OLD vs SLOVENIA comparison ===")
print(f"Samples compared: {len(rows)}")
print()
print(f"Same top species: {same_species}/{len(rows)}")
print(f"Same top genus:   {same_genus}/{len(rows)}")
print(f"Same top family:  {same_family}/{len(rows)}")
print()
print("Created:")
print(output_file)
