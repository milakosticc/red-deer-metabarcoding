import csv

qc_file = "taxonomy/qc/LIB1_taxonomy_QC.tsv"
abundance_file = "taxonomy/input/LIB1_final_ASV_table.csv"

out_table = "taxonomy/final/LIB1_target_plant_ASV_table.csv"
out_ids = "taxonomy/final/LIB1_target_plant_ASV_ids.txt"
out_taxonomy = "taxonomy/final/LIB1_target_plant_abundance_taxonomy.tsv"
summary_file = "taxonomy/final/LIB1_target_filter_summary.txt"
regional_qc = "taxonomy/qc/LIB1_species_for_regional_QC.tsv"

# Read QC table
target_ids = []
qc_rows = {}

with open(qc_file, encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    qc_fields = reader.fieldnames

    for row in reader:
        qc_rows[row["ASV_ID"]] = row

        if row["QC_status"] == "TARGET_PLANT":
            target_ids.append(row["ASV_ID"])

target_set = set(target_ids)

# Save target IDs
with open(out_ids, "w") as f:
    for asv in target_ids:
        f.write(asv + "\n")

# Filter abundance table
with open(abundance_file, encoding="utf-8", newline="") as f:
    reader = csv.reader(f)

    header = next(reader)
    all_asvs = header[1:]

    keep_indices = [
        i for i, asv in enumerate(all_asvs)
        if asv in target_set
    ]

    kept_asvs = [all_asvs[i] for i in keep_indices]

    total_target_reads = 0

    with open(out_table, "w", newline="") as out:
        writer = csv.writer(out)

        writer.writerow([""] + kept_asvs)

        for row in reader:
            sample = row[0]
            counts = [int(float(row[i + 1])) for i in keep_indices]

            total_target_reads += sum(counts)

            writer.writerow([sample] + counts)

# Write target abundance + taxonomy table
target_rows = [
    qc_rows[x]
    for x in target_ids
]

target_rows.sort(
    key=lambda x: int(x["Total_reads"]),
    reverse=True
)

with open(out_taxonomy, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=qc_fields,
        delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(target_rows)

# Regional QC: only confident species assignments
species_rows = [
    x for x in target_rows
    if x.get("Species", "").strip()
]

regional_fields = [
    "ASV_ID",
    "Total_reads",
    "Prevalence_samples",
    "Family",
    "Genus",
    "Species",
    "Species_confidence"
]

with open(regional_qc, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=regional_fields,
        delimiter="\t",
        extrasaction="ignore"
    )
    writer.writeheader()
    writer.writerows(species_rows)

# Summary
with open(summary_file, "w") as f:
    f.write("=== LIB1 target-taxa filtering ===\n")
    f.write(f"Input ASVs: {len(qc_rows)}\n")
    f.write(f"Target plant ASVs retained: {len(target_ids)}\n")
    f.write(f"Target plant reads retained: {total_target_reads}\n")
    f.write(f"ASV retention: {100*len(target_ids)/len(qc_rows):.2f}%\n")
    f.write(f"Read retention: {100*total_target_reads/473081:.2f}%\n")

print("Target filtering completed.")
print("Target ASVs:", len(target_ids))
print("Target reads:", total_target_reads)
