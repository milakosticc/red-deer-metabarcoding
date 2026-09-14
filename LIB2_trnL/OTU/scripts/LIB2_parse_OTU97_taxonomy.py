import csv
from collections import Counter

input_file = "taxonomy/otu_97/LIB2_OTU97_SINTAX_0.8.tsv"
output_file = "taxonomy/otu_97/LIB2_OTU97_taxonomy.tsv"

ranks = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species"
}

rank_order = [
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
]

rows = []

with open(input_file, encoding="utf-8") as f:
    for line in f:
        if not line.strip():
            continue

        parts = line.rstrip("\n").split("\t")

        query = parts[0]

        if "|centroid=" in query:
            otu_id, centroid = query.split("|centroid=", 1)
        else:
            otu_id = query
            centroid = ""

        # Column 4 contains taxonomy accepted at the SINTAX cutoff
        accepted_taxonomy = parts[3] if len(parts) >= 4 else ""

        taxonomy = {rank: "" for rank in rank_order}

        if accepted_taxonomy:
            for item in accepted_taxonomy.split(","):
                item = item.strip()

                if ":" not in item:
                    continue

                prefix, name = item.split(":", 1)

                if prefix in ranks:
                    taxonomy[ranks[prefix]] = name.replace("_", " ")

        deepest = "Unassigned"

        for rank in rank_order:
            if taxonomy[rank]:
                deepest = rank

        rows.append({
            "OTU_ID": otu_id,
            "Centroid_ASV": centroid,
            **taxonomy,
            "Deepest_Rank": deepest
        })

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=[
            "OTU_ID",
            "Centroid_ASV",
            *rank_order,
            "Deepest_Rank"
        ],
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(rows)

print("=== LIB2 OTU97 TAXONOMY SUMMARY ===")
print("Total OTUs:", len(rows))

for rank in rank_order:
    n = sum(bool(row[rank]) for row in rows)
    print(f"At least {rank}: {n}")

deepest_counts = Counter(row["Deepest_Rank"] for row in rows)

print("\nDeepest assigned rank:")
for rank in reversed(rank_order):
    print(f"{rank}: {deepest_counts.get(rank, 0)}")

print("Unassigned:", deepest_counts.get("Unassigned", 0))
print("\nOutput:", output_file)
