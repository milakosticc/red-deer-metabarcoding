from collections import Counter

inp = "taxonomy/otu_97/LIB1_OTU97_BOLD_sintax_0.8.tsv"
out = "taxonomy/otu_97/LIB1_OTU97_taxonomy.tsv"

ranks = [
    ("k", "Kingdom"),
    ("p", "Phylum"),
    ("c", "Class"),
    ("o", "Order"),
    ("f", "Family"),
    ("g", "Genus"),
    ("s", "Species")
]

rows = []
deepest_counts = Counter()

with open(inp, encoding="utf-8") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")

        query = parts[0]
        otu = query.split("|")[0]

        confident = parts[3] if len(parts) >= 4 else ""

        tax = {}

        if confident:
            for x in confident.split(","):
                if ":" in x:
                    prefix, value = x.split(":", 1)
                    tax[prefix] = value

        deepest = "Unassigned"

        for prefix, name in ranks:
            if tax.get(prefix):
                deepest = name

        deepest_counts[deepest] += 1

        rows.append([
            otu,
            tax.get("k", ""),
            tax.get("p", ""),
            tax.get("c", ""),
            tax.get("o", ""),
            tax.get("f", ""),
            tax.get("g", ""),
            tax.get("s", ""),
            deepest
        ])

with open(out, "w", encoding="utf-8") as f:
    f.write(
        "OTU_ID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tDeepest_confident_rank\n"
    )

    for row in rows:
        f.write("\t".join(row) + "\n")

print("=== LIB1 OTU97 TAXONOMY SUMMARY ===")
print("Total OTUs:", len(rows))

for i, (_, name) in enumerate(ranks, start=1):
    n = sum(1 for r in rows if r[i])
    print(f"At least {name}: {n}")

print("\n=== Deepest confident rank ===")

for name in [
    "Species", "Genus", "Family",
    "Order", "Class", "Phylum",
    "Kingdom", "Unassigned"
]:
    print(f"{name}: {deepest_counts[name]}")

print("\nOutput:", out)
