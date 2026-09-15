#!/usr/bin/env python3

from pathlib import Path

input_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_OTU97_Slovenia_BOLD_sintax_0.8.tsv"
)

output_file = Path(
    "taxonomy/otu_97/slovenia/LIB1_OTU97_Slovenia_taxonomy.tsv"
)

ranks = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species"
}

rank_order = ["k", "p", "c", "o", "f", "g", "s"]

with open(input_file) as infile, open(output_file, "w") as outfile:

    header = [
        "OTU",
        "Centroid_ASV",
        "Kingdom",
        "Phylum",
        "Class",
        "Order",
        "Family",
        "Genus",
        "Species",
        "Deepest_level"
    ]

    outfile.write("\t".join(header) + "\n")

    for line in infile:
        parts = line.rstrip("\n").split("\t")

        query = parts[0]
        accepted_taxonomy = parts[3] if len(parts) > 3 else ""

        # OTU1|centroid=ASV398
        query_parts = query.split("|")
        otu = query_parts[0]

        centroid = ""
        for x in query_parts[1:]:
            if x.startswith("centroid="):
                centroid = x.split("=", 1)[1]

        taxonomy = {rank: "" for rank in rank_order}

        for item in accepted_taxonomy.split(","):
            if ":" not in item:
                continue

            rank, value = item.split(":", 1)

            if rank in taxonomy:
                taxonomy[rank] = value.replace("_", " ")

        deepest = ""

        for rank in rank_order:
            if taxonomy[rank]:
                deepest = ranks[rank]

        row = [
            otu,
            centroid,
            taxonomy["k"],
            taxonomy["p"],
            taxonomy["c"],
            taxonomy["o"],
            taxonomy["f"],
            taxonomy["g"],
            taxonomy["s"],
            deepest
        ]

        outfile.write("\t".join(row) + "\n")

print("Final Slovenia taxonomy table created:")
print(output_file)
