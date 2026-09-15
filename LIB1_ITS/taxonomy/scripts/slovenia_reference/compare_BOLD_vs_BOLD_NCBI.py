#!/usr/bin/env python3

from pathlib import Path

bold_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_OTU97_Slovenia_BOLD_sintax_0.8.tsv"
)

combined_file = Path(
    "taxonomy/otu_97/slovenia_bold_ncbi/"
    "LIB1_OTU97_Slovenia_BOLD_NCBI_sintax_0.8.tsv"
)


def read_taxonomy(path):
    result = {}

    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            otu = parts[0].split("|")[0]

            accepted = parts[3] if len(parts) > 3 else ""

            result[otu] = accepted

    return result


def deepest_rank(tax):
    ranks = [
        ("s:", "Species"),
        ("g:", "Genus"),
        ("f:", "Family"),
        ("o:", "Order"),
        ("c:", "Class"),
        ("p:", "Phylum"),
        ("k:", "Kingdom")
    ]

    for prefix, name in ranks:
        if tax.startswith(prefix) or f",{prefix}" in tax:
            return name

    return "Unclassified"


bold = read_taxonomy(bold_file)
combined = read_taxonomy(combined_file)

same_taxonomy = 0
changed_taxonomy = 0
same_deepest = 0
changed_deepest = 0

changes = []

for otu in sorted(bold):

    old = bold[otu]
    new = combined.get(otu, "")

    old_level = deepest_rank(old)
    new_level = deepest_rank(new)

    if old == new:
        same_taxonomy += 1
    else:
        changed_taxonomy += 1
        changes.append(
            (otu, old_level, new_level, old, new)
        )

    if old_level == new_level:
        same_deepest += 1
    else:
        changed_deepest += 1


print("=== BOLD-only vs BOLD+NCBI ===")
print(f"OTUs compared: {len(bold)}")
print(f"Identical accepted taxonomy: {same_taxonomy}")
print(f"Changed accepted taxonomy:   {changed_taxonomy}")
print()
print(f"Same deepest rank:            {same_deepest}")
print(f"Changed deepest rank:         {changed_deepest}")

print()
print("=== CHANGED OTUs ===")

for otu, old_level, new_level, old, new in changes:
    print()
    print(otu)
    print(f"  BOLD-only    [{old_level}]: {old}")
    print(f"  BOLD+NCBI    [{new_level}]: {new}")
