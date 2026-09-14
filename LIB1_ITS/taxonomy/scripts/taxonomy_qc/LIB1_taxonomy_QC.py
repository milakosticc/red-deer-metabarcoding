import csv

abundance_tax = "taxonomy/final/LIB1_ASV_abundance_taxonomy.tsv"
euk_file = "taxonomy/sintax/LIB1_EUKARYOME_sintax_0.8.tsv"

out_file = "taxonomy/qc/LIB1_taxonomy_QC.tsv"
summary_file = "taxonomy/qc/LIB1_taxonomy_QC_summary.txt"
review_file = "taxonomy/qc/LIB1_ASVs_for_review.tsv"

# -----------------------------
# Read EUKARYOME classification
# -----------------------------
euk = {}

with open(euk_file, encoding="utf-8") as f:
    for line in f:
        parts = line.rstrip("\n").split("\t")

        asv = parts[0]

        # VSEARCH output has 4 columns
        confident = parts[3].strip() if len(parts) >= 4 else ""

        if confident.startswith("d:Viridiplantae"):
            status = "TARGET_PLANT"

        elif confident == "":
            status = "UNASSIGNED"

        else:
            status = "NON_TARGET"

        euk[asv] = {
            "EUKARYOME_taxonomy": confident,
            "QC_status": status
        }

# -----------------------------
# Merge with abundance + BOLD taxonomy
# -----------------------------
records = []

with open(abundance_tax, encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")

    original_fields = reader.fieldnames

    for row in reader:
        asv = row["ASV_ID"]

        row["QC_status"] = euk.get(asv, {}).get(
            "QC_status", "UNASSIGNED"
        )

        row["EUKARYOME_taxonomy"] = euk.get(asv, {}).get(
            "EUKARYOME_taxonomy", ""
        )

        records.append(row)

# -----------------------------
# Write full QC table
# -----------------------------
out_fields = (
    ["ASV_ID", "Total_reads", "Prevalence_samples",
     "QC_status", "EUKARYOME_taxonomy"]
    +
    [x for x in original_fields
     if x not in ["ASV_ID", "Total_reads", "Prevalence_samples"]]
)

with open(out_file, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=out_fields,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(records)

# -----------------------------
# Summary
# -----------------------------
statuses = ["TARGET_PLANT", "NON_TARGET", "UNASSIGNED"]

total_asvs = len(records)
total_reads = sum(int(x["Total_reads"]) for x in records)

with open(summary_file, "w", encoding="utf-8") as f:
    f.write("=== LIB1 taxonomy QC summary ===\n")
    f.write(f"Total ASVs: {total_asvs}\n")
    f.write(f"Total reads: {total_reads}\n\n")

    for status in statuses:

        subset = [
            x for x in records
            if x["QC_status"] == status
        ]

        n_asv = len(subset)
        n_reads = sum(int(x["Total_reads"]) for x in subset)

        asv_pct = 100 * n_asv / total_asvs
        read_pct = 100 * n_reads / total_reads if total_reads else 0

        f.write(
            f"{status}: "
            f"{n_asv} ASVs ({asv_pct:.2f}%), "
            f"{n_reads} reads ({read_pct:.2f}%)\n"
        )

# -----------------------------
# Important ASVs requiring review
# -----------------------------
review = [
    x for x in records
    if x["QC_status"] != "TARGET_PLANT"
]

review.sort(
    key=lambda x: int(x["Total_reads"]),
    reverse=True
)

review_fields = [
    "ASV_ID",
    "Total_reads",
    "Prevalence_samples",
    "QC_status",
    "EUKARYOME_taxonomy",
    "Family",
    "Genus",
    "Species",
    "Deepest_confident_rank"
]

with open(review_file, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=review_fields,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(review)

print("Created:", out_file)
print("Created:", summary_file)
print("Created:", review_file)
