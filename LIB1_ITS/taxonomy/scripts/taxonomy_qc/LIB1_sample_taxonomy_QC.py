import csv

abundance_file = "taxonomy/input/LIB1_final_ASV_table.csv"
qc_file = "taxonomy/qc/LIB1_taxonomy_QC.tsv"

out_file = "taxonomy/qc/LIB1_sample_taxonomy_QC.tsv"

# Read ASV -> QC status
status = {}

with open(qc_file, encoding="utf-8", newline="") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        status[row["ASV_ID"]] = row["QC_status"]

# Read abundance table
with open(abundance_file, encoding="utf-8", newline="") as f:
    reader = csv.reader(f)

    header = next(reader)
    asvs = header[1:]

    results = []

    for row in reader:
        sample = row[0]
        counts = list(map(int, row[1:]))

        totals = {
            "TARGET_PLANT": 0,
            "NON_TARGET": 0,
            "UNASSIGNED": 0
        }

        for asv, count in zip(asvs, counts):
            s = status.get(asv, "UNASSIGNED")
            totals[s] += count

        total_reads = sum(counts)

        results.append({
            "Sample": sample,
            "Total_reads": total_reads,
            "Plant_reads": totals["TARGET_PLANT"],
            "Plant_percent": 100 * totals["TARGET_PLANT"] / total_reads if total_reads else 0,
            "Non_target_reads": totals["NON_TARGET"],
            "Non_target_percent": 100 * totals["NON_TARGET"] / total_reads if total_reads else 0,
            "Unassigned_reads": totals["UNASSIGNED"],
            "Unassigned_percent": 100 * totals["UNASSIGNED"] / total_reads if total_reads else 0
        })

fields = [
    "Sample",
    "Total_reads",
    "Plant_reads",
    "Plant_percent",
    "Non_target_reads",
    "Non_target_percent",
    "Unassigned_reads",
    "Unassigned_percent"
]

with open(out_file, "w", encoding="utf-8", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    writer.writeheader()

    for r in results:
        r["Plant_percent"] = f'{r["Plant_percent"]:.2f}'
        r["Non_target_percent"] = f'{r["Non_target_percent"]:.2f}'
        r["Unassigned_percent"] = f'{r["Unassigned_percent"]:.2f}'
        writer.writerow(r)

print("Created:", out_file)
