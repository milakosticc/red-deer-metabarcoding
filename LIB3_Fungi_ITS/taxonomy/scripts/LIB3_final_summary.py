import csv

tax_file = "taxonomy/final/LIB3_target_fungi_abundance_taxonomy.tsv"
seq_tracking = "dada2_out/final/LIB3_final_sequence_tracking.csv"
asv_tracking = "dada2_out/final/LIB3_final_ASV_tracking.csv"

out_seq = "taxonomy/final/LIB3_final_sequence_tracking.tsv"
out_asv = "taxonomy/final/LIB3_final_ASV_tracking.tsv"
out_summary = "taxonomy/final/LIB3_final_taxonomy_summary.txt"

# -------------------------
# Final fungal taxonomy
# -------------------------
with open(tax_file, newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

total_asvs = len(rows)
total_reads = sum(int(r["Total_reads"]) for r in rows)

family_rows = [r for r in rows if r["Family"]]
genus_rows = [r for r in rows if r["Genus"]]
species_rows = [r for r in rows if r["Species"]]

families = {r["Family"] for r in family_rows}
genera = {r["Genus"] for r in genus_rows}
species = {r["Species"] for r in species_rows}

# -------------------------
# Sequence tracking
# -------------------------
with open(seq_tracking, newline="", encoding="utf-8") as f:
    seq_rows = list(csv.DictReader(f))

with open(out_seq, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Step", "Reads", "Status"])

    for r in seq_rows:
        w.writerow([
            r["Step"],
            r["Reads"],
            r["Status"]
        ])

    w.writerow([
        "Target-fungi-filtered",
        total_reads,
        "Completed"
    ])

# -------------------------
# ASV tracking
# -------------------------
with open(asv_tracking, newline="", encoding="utf-8") as f:
    asv_rows = list(csv.DictReader(f))

with open(out_asv, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Step", "ASVs"])

    for r in asv_rows:
        w.writerow([
            r["Step"],
            r["ASVs"]
        ])

    w.writerow([
        "Target-fungi-filtered",
        total_asvs
    ])

# -------------------------
# Summary
# -------------------------
summary = [
    f"Final fungal ASVs: {total_asvs}",
    f"Final fungal reads: {total_reads}",
    f"ASVs with family assignment: {len(family_rows)}",
    f"Unique families: {len(families)}",
    f"ASVs with genus assignment: {len(genus_rows)}",
    f"Unique genera: {len(genera)}",
    f"ASVs with species assignment: {len(species_rows)}",
    f"Unique species-level taxa: {len(species)}"
]

with open(out_summary, "w", encoding="utf-8") as f:
    f.write("\n".join(summary) + "\n")

print("=== LIB3 FINAL TAXONOMY SUMMARY ===")
for x in summary:
    print(x)

print("\n=== FINAL SEQUENCE TRACKING ===")
with open(out_seq) as f:
    print(f.read())

print("=== FINAL ASV TRACKING ===")
with open(out_asv) as f:
    print(f.read())
