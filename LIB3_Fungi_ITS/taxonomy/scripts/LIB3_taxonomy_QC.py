import csv
from collections import Counter

euk_file = "taxonomy/sintax/LIB3_EUKARYOME_sintax_0.8.tsv"
unite_file = "taxonomy/final/LIB3_taxonomy_UNITE_SINTAX_0.8.tsv"
ab_file = "taxonomy/input/LIB3_final_ASV_table.csv"

out_qc = "taxonomy/qc/LIB3_taxonomy_QC.tsv"
out_review = "taxonomy/qc/LIB3_ASVs_for_review.tsv"

# UNITE taxonomy
unite = {}
with open(unite_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        unite[row["ASV_ID"]] = row

asv_ids = set(unite)

# EUKARYOME broad classification
euk = {}
with open(euk_file, encoding="utf-8") as f:
    for line in f:
        p = line.rstrip("\n").split("\t")
        asv = p[0]
        confident = p[3] if len(p) >= 4 else ""

        if not confident:
            status = "UNASSIGNED"
        elif confident.startswith("d:Fungi"):
            status = "TARGET_FUNGI"
        else:
            status = "NON_TARGET"

        euk[asv] = (status, confident)

# Abundance table
with open(ab_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [c for c in headers if c in asv_ids]

if not asv_cols:
    raise ValueError("ASV columns not found in abundance table.")

records = []

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

    status, assignment = euk.get(asv, ("UNASSIGNED", ""))
    u = unite[asv]

    records.append({
        "ASV_ID": asv,
        "Total_reads": int(total),
        "Prevalence_samples": prevalence,
        "QC_status": status,
        "EUKARYOME_assignment": assignment,
        "Phylum": u.get("Phylum", ""),
        "Class": u.get("Class", ""),
        "Order": u.get("Order", ""),
        "Family": u.get("Family", ""),
        "Genus": u.get("Genus", ""),
        "Species": u.get("Species", ""),
        "Deepest_confident_rank": u.get("Deepest_confident_rank", "")
    })

records.sort(key=lambda x: x["Total_reads"], reverse=True)

fields = list(records[0].keys())

with open(out_qc, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows(records)

with open(out_review, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
    w.writeheader()
    w.writerows([r for r in records if r["QC_status"] != "TARGET_FUNGI"])

counts = Counter(r["QC_status"] for r in records)
reads = Counter()

for r in records:
    reads[r["QC_status"]] += r["Total_reads"]

total_reads = sum(r["Total_reads"] for r in records)

print("=== LIB3 BROAD TAXONOMY QC ===")
print("Total ASVs:", len(records))
print("Total reads:", total_reads)
print()

for status in ["TARGET_FUNGI", "NON_TARGET", "UNASSIGNED"]:
    n = counts[status]
    rr = reads[status]

    print(
        f"{status}: {n} ASVs ({n/len(records)*100:.2f}%), "
        f"{rr} reads ({rr/total_reads*100:.2f}%)"
    )
