import csv
from collections import defaultdict

table_file = "taxonomy/final/LIB3_target_fungi_ASV_table.csv"
tax_file = "taxonomy/final/LIB3_target_fungi_abundance_taxonomy.tsv"
out_file = "taxonomy/final/LIB3_top5_species_genus_family_per_sample.tsv"

# -------------------------
# Read taxonomy
# -------------------------
taxonomy = {}

with open(tax_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        taxonomy[row["ASV_ID"]] = {
            "Species": row["Species"],
            "Genus": row["Genus"],
            "Family": row["Family"]
        }

asv_ids = set(taxonomy)

# -------------------------
# Read final fungal table
# -------------------------
with open(table_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [x for x in headers if x in asv_ids]
metadata_cols = [x for x in headers if x not in asv_ids]

if not metadata_cols:
    raise ValueError("Could not identify sample column.")

sample_col = metadata_cols[0]

print("Sample column:", sample_col)
print("Samples:", len(rows))
print("Fungal ASV columns:", len(asv_cols))

# -------------------------
# Aggregate taxa per sample
# -------------------------
output = []

for row in rows:

    sample = row[sample_col]

    asv_counts = {}

    for asv in asv_cols:
        try:
            n = int(float(row.get(asv, 0) or 0))
        except ValueError:
            n = 0

        asv_counts[asv] = n

    total_reads = sum(asv_counts.values())

    record = {
        "Sample": sample,
        "Total_fungal_reads": total_reads
    }

    for rank in ["Species", "Genus", "Family"]:

        totals = defaultdict(int)

        for asv, n in asv_counts.items():

            if n <= 0:
                continue

            taxon = taxonomy[asv][rank].strip()

            if taxon:
                totals[taxon] += n

        ranked = sorted(
            totals.items(),
            key=lambda x: x[1],
            reverse=True
        )[:5]

        for i in range(1, 6):

            if i <= len(ranked):
                taxon, reads = ranked[i - 1]

                pct = (
                    reads / total_reads * 100
                    if total_reads > 0 else 0
                )

                record[f"Top_{i}_{rank}"] = (
                    f"{taxon} | {reads} reads | {pct:.2f}%"
                )

            else:
                record[f"Top_{i}_{rank}"] = ""

    output.append(record)

# -------------------------
# Write output
# -------------------------
fields = ["Sample", "Total_fungal_reads"]

for rank in ["Species", "Genus", "Family"]:
    for i in range(1, 6):
        fields.append(f"Top_{i}_{rank}")

with open(out_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=fields,
        delimiter="\t"
    )

    writer.writeheader()
    writer.writerows(output)

print("\nCreated:", out_file)

print("\n=== FIRST 3 SAMPLES ===")

for r in output[:3]:

    print("\nSample:", r["Sample"])
    print("Total fungal reads:", r["Total_fungal_reads"])

    for rank in ["Species", "Genus", "Family"]:

        print(rank + ":")

        for i in range(1, 6):

            value = r[f"Top_{i}_{rank}"]

            if value:
                print(" ", i, value)
