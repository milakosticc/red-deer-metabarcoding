import csv

abundance_file = "taxonomy/input/LIB2_final_ASV_table.csv"
taxonomy_file = "taxonomy/final/LIB2_taxonomy_SINTAX_0.8.tsv"
output_file = "taxonomy/final/LIB2_top5_species_genus_family_per_sample.tsv"

# Read taxonomy
taxonomy = {}

with open(taxonomy_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        taxonomy[row["ASV_ID"]] = row

asv_ids = set(taxonomy)

# Read abundance table
with open(abundance_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [c for c in headers if c in asv_ids]

if not asv_cols:
    raise ValueError("ASV columns not found in abundance table.")

# Detect sample column
non_asv_cols = [c for c in headers if c not in asv_cols]

sample_col = None
for c in non_asv_cols:
    values = [str(r.get(c, "")) for r in rows]
    if any(v.startswith("LME") for v in values):
        sample_col = c
        break

if sample_col is None:
    sample_col = non_asv_cols[0]

def valid_taxon(x):
    return x not in ("", "NA", "nan", "None", None)

def clean_name(x):
    return str(x).replace("_", " ")

def top5_for_rank(row, rank, total_reads):
    counts = {}

    for asv in asv_cols:
        tax = taxonomy.get(asv, {}).get(rank, "")

        if not valid_taxon(tax):
            continue

        try:
            reads = int(float(row.get(asv, 0) or 0))
        except ValueError:
            reads = 0

        if reads <= 0:
            continue

        counts[tax] = counts.get(tax, 0) + reads

    ranked = sorted(
        counts.items(),
        key=lambda x: x[1],
        reverse=True
    )[:5]

    result = []

    for taxon, reads in ranked:
        pct = (reads / total_reads * 100) if total_reads > 0 else 0
        result.append(
            f"{clean_name(taxon)} ({reads:,}; {pct:.2f}%)"
        )

    while len(result) < 5:
        result.append("—")

    return result


output_rows = []

for row in rows:
    sample = row[sample_col]

    total_reads = 0

    for asv in asv_cols:
        try:
            total_reads += int(float(row.get(asv, 0) or 0))
        except ValueError:
            pass

    species = top5_for_rank(row, "Species", total_reads)
    genus   = top5_for_rank(row, "Genus", total_reads)
    family  = top5_for_rank(row, "Family", total_reads)

    output_rows.append({
        "Sample": sample,

        "Top_1_Species": species[0],
        "Top_2_Species": species[1],
        "Top_3_Species": species[2],
        "Top_4_Species": species[3],
        "Top_5_Species": species[4],

        "Top_1_Genus": genus[0],
        "Top_2_Genus": genus[1],
        "Top_3_Genus": genus[2],
        "Top_4_Genus": genus[3],
        "Top_5_Genus": genus[4],

        "Top_1_Family": family[0],
        "Top_2_Family": family[1],
        "Top_3_Family": family[2],
        "Top_4_Family": family[3],
        "Top_5_Family": family[4]
    })

fieldnames = [
    "Sample",
    "Top_1_Species",
    "Top_2_Species",
    "Top_3_Species",
    "Top_4_Species",
    "Top_5_Species",
    "Top_1_Genus",
    "Top_2_Genus",
    "Top_3_Genus",
    "Top_4_Genus",
    "Top_5_Genus",
    "Top_1_Family",
    "Top_2_Family",
    "Top_3_Family",
    "Top_4_Family",
    "Top_5_Family"
]

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    writer.writeheader()
    writer.writerows(output_rows)

print("Created:", output_file)
print("Samples:", len(output_rows))
