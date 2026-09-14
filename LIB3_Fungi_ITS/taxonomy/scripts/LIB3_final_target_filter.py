import csv

qc_file = "taxonomy/qc/LIB3_taxonomy_QC.tsv"
table_file = "taxonomy/input/LIB3_final_ASV_table.csv"
fasta_file = "taxonomy/input/LIB3_final_ASVs.fasta"

out_table = "taxonomy/final/LIB3_target_fungi_ASV_table.csv"
out_fasta = "taxonomy/final/LIB3_target_fungi_ASVs.fasta"
out_tax = "taxonomy/final/LIB3_target_fungi_abundance_taxonomy.tsv"
out_summary = "taxonomy/final/LIB3_target_filter_summary.txt"

# -------------------------
# Read QC / taxonomy
# -------------------------
qc = {}

with open(qc_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        qc[row["ASV_ID"]] = row

target = {
    asv for asv, row in qc.items()
    if row["QC_status"] == "TARGET_FUNGI"
}

# -------------------------
# Read abundance table
# -------------------------
with open(table_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    fields = reader.fieldnames

asv_cols = [x for x in fields if x in qc]
metadata_cols = [x for x in fields if x not in qc]

target_cols = [x for x in asv_cols if x in target]

# Write filtered abundance table
with open(out_table, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=metadata_cols + target_cols
    )
    writer.writeheader()

    for row in rows:
        writer.writerow({
            k: row[k]
            for k in metadata_cols + target_cols
        })

# -------------------------
# Calculate total abundance
# -------------------------
totals = {}

for asv in target_cols:
    total = 0

    for row in rows:
        try:
            total += int(float(row.get(asv, 0) or 0))
        except ValueError:
            pass

    totals[asv] = total

# -------------------------
# Write abundance + taxonomy
# -------------------------
tax_fields = [
    "ASV_ID",
    "Total_reads",
    "Prevalence_samples",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species",
    "Deepest_confident_rank"
]

with open(out_tax, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=tax_fields,
        delimiter="\t"
    )
    writer.writeheader()

    for asv in sorted(
        target_cols,
        key=lambda x: totals[x],
        reverse=True
    ):
        r = qc[asv]

        writer.writerow({
            "ASV_ID": asv,
            "Total_reads": totals[asv],
            "Prevalence_samples": r["Prevalence_samples"],
            "Phylum": r["Phylum"],
            "Class": r["Class"],
            "Order": r["Order"],
            "Family": r["Family"],
            "Genus": r["Genus"],
            "Species": r["Species"],
            "Deepest_confident_rank":
                r["Deepest_confident_rank"]
        })

# -------------------------
# Filter FASTA
# -------------------------
written = 0
keep = False

with open(fasta_file, encoding="utf-8") as inp, \
     open(out_fasta, "w", encoding="utf-8") as out:

    for line in inp:
        if line.startswith(">"):
            asv = line[1:].strip().split()[0]
            keep = asv in target

            if keep:
                written += 1

        if keep:
            out.write(line)

# -------------------------
# Summary
# -------------------------
target_reads = sum(totals.values())
all_reads = sum(
    int(r["Total_reads"])
    for r in qc.values()
)

removed_reads = all_reads - target_reads

with open(out_summary, "w") as f:
    f.write("LIB3 FINAL TARGET FILTER\n")
    f.write(f"Input ASVs: {len(qc)}\n")
    f.write(f"Target fungal ASVs: {len(target)}\n")
    f.write(f"Removed ASVs: {len(qc)-len(target)}\n")
    f.write(f"Input reads: {all_reads}\n")
    f.write(f"Target fungal reads: {target_reads}\n")
    f.write(f"Removed reads: {removed_reads}\n")
    f.write(
        f"Reads retained: "
        f"{target_reads/all_reads*100:.2f}%\n"
    )

print("=== LIB3 FINAL TARGET FILTER ===")
print("Input ASVs:", len(qc))
print("Target fungal ASVs:", len(target))
print("FASTA sequences written:", written)
print("Input reads:", all_reads)
print("Target fungal reads:", target_reads)
print("Removed reads:", removed_reads)
print(
    "Reads retained:",
    f"{target_reads/all_reads*100:.2f}%"
)
print()
print("Final table:", out_table)
print("Final FASTA:", out_fasta)
print("Taxonomy:", out_tax)
