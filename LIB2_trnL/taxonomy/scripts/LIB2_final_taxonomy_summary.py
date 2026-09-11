import csv

tax_file = "taxonomy/final/LIB2_taxonomy_SINTAX_0.8.tsv"
ab_file  = "taxonomy/input/LIB2_final_ASV_table.csv"

# -----------------------------
# Read taxonomy
# -----------------------------
taxonomy = {}

with open(tax_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        taxonomy[row["ASV_ID"]] = row

asv_ids = set(taxonomy)

# -----------------------------
# Read abundance table
# -----------------------------
with open(ab_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [x for x in headers if x in asv_ids]

records = []

# Normal format: samples = rows, ASVs = columns
if asv_cols:
    for asv in asv_cols:
        total = 0
        prevalence = 0

        for row in rows:
            try:
                value = float(row.get(asv, 0) or 0)
            except ValueError:
                value = 0

            total += value

            if value > 0:
                prevalence += 1

        r = {
            "ASV_ID": asv,
            "Total_reads": int(total),
            "Prevalence_samples": prevalence
        }

        r.update(taxonomy[asv])
        records.append(r)

else:
    raise ValueError("ASV columns could not be found in abundance table.")

records.sort(key=lambda x: x["Total_reads"], reverse=True)

# -----------------------------
# Save abundance + taxonomy
# -----------------------------
columns = [
    "ASV_ID",
    "Total_reads",
    "Prevalence_samples",
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Deepest_confident_rank"
]

with open(
    "taxonomy/final/LIB2_ASV_abundance_taxonomy.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=columns,
        delimiter="\t",
        extrasaction="ignore"
    )

    writer.writeheader()
    writer.writerows(records)

# Top 30
with open(
    "taxonomy/qc/LIB2_top30_ASVs.tsv",
    "w",
    newline="",
    encoding="utf-8"
) as f:

    writer = csv.DictWriter(
        f,
        fieldnames=columns,
        delimiter="\t",
        extrasaction="ignore"
    )

    writer.writeheader()
    writer.writerows(records[:30])


def valid(x):
    return x not in ("", "NA", "nan", "None", None)


def assigned(rank):
    return sum(1 for r in records if valid(r.get(rank)))


def unique(rank):
    return len({
        r.get(rank)
        for r in records
        if valid(r.get(rank))
    })


print("=== LIB2 FINAL TAXONOMY SUMMARY ===")
print("Total ASVs:", len(records))
print("Total reads:", sum(r["Total_reads"] for r in records))

print()
print("ASVs with family assignment:", assigned("Family"))
print("Unique families:", unique("Family"))

print()
print("ASVs with genus assignment:", assigned("Genus"))
print("Unique genera:", unique("Genus"))

print()
print("ASVs with species assignment:", assigned("Species"))
print("Unique species-level taxa:", unique("Species"))

print("\n=== TOP 15 ASVs ===")

print(
    "ASV_ID\tReads\tPrev\tFamily\tGenus\tSpecies\tDeepest"
)

for r in records[:15]:
    print(
        f'{r["ASV_ID"]}\t'
        f'{r["Total_reads"]}\t'
        f'{r["Prevalence_samples"]}\t'
        f'{r.get("Family","")}\t'
        f'{r.get("Genus","")}\t'
        f'{r.get("Species","")}\t'
        f'{r.get("Deepest_confident_rank","")}'
    )
