import csv
from collections import defaultdict

table_file = "taxonomy/input/LIB2_final_ASV_table.csv"
tax_file = "taxonomy/final/LIB2_ASV_abundance_taxonomy.tsv"
outdir = "taxonomy/ecology_input"

# Read taxonomy
taxonomy = {}

with open(tax_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")

    for row in reader:
        taxonomy[row["ASV_ID"]] = {
            "Family": row.get("Family", "").strip(),
            "Genus": row.get("Genus", "").strip(),
            "Species": row.get("Species", "").strip()
        }

asv_ids = set(taxonomy)

# Read abundance table
with open(table_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [c for c in headers if c in asv_ids]
metadata_cols = [c for c in headers if c not in asv_ids]

if not metadata_cols:
    raise ValueError("No sample/metadata column detected.")

sample_col = metadata_cols[0]

# Sample totals
sample_totals = {}

for row in rows:
    sample = row[sample_col]

    total = 0

    for asv in asv_cols:
        try:
            total += int(float(row.get(asv, 0) or 0))
        except ValueError:
            pass

    sample_totals[sample] = total

with open(f"{outdir}/LIB2_sample_read_totals.tsv",
          "w", newline="", encoding="utf-8") as f:

    w = csv.writer(f, delimiter="\t")
    w.writerow(["Sample", "Total_LIB2_reads"])

    for sample, total in sample_totals.items():
        w.writerow([sample, total])

# Aggregate at each taxonomic rank
for rank in ["Family", "Genus", "Species"]:

    sample_taxa = {}
    all_taxa = set()

    for row in rows:
        sample = row[sample_col]
        counts = defaultdict(int)

        for asv in asv_cols:
            try:
                n = int(float(row.get(asv, 0) or 0))
            except ValueError:
                n = 0

            if n == 0:
                continue

            taxon = taxonomy[asv][rank]

            if not taxon:
                taxon = f"Unclassified_at_{rank.lower()}"

            counts[taxon] += n
            all_taxa.add(taxon)

        sample_taxa[sample] = counts

    taxa = sorted(all_taxa)

    # Count table
    with open(
        f"{outdir}/LIB2_{rank.lower()}_abundance_counts.tsv",
        "w", newline="", encoding="utf-8"
    ) as f:

        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample"] + taxa)

        for sample in sample_totals:
            w.writerow(
                [sample] +
                [sample_taxa[sample].get(t, 0) for t in taxa]
            )

    # Relative abundance table
    with open(
        f"{outdir}/LIB2_{rank.lower()}_relative_abundance.tsv",
        "w", newline="", encoding="utf-8"
    ) as f:

        w = csv.writer(f, delimiter="\t")
        w.writerow(["Sample"] + taxa)

        for sample in sample_totals:
            total = sample_totals[sample]

            values = []

            for taxon in taxa:
                n = sample_taxa[sample].get(taxon, 0)

                pct = (
                    n / total * 100
                    if total > 0 else 0
                )

                values.append(f"{pct:.6f}")

            w.writerow([sample] + values)

    unclassified = f"Unclassified_at_{rank.lower()}"

    classified_reads = 0
    unclassified_reads = 0

    for sample in sample_taxa:
        for taxon, n in sample_taxa[sample].items():

            if taxon == unclassified:
                unclassified_reads += n
            else:
                classified_reads += n

    print(f"\n=== {rank.upper()} ===")
    print("Taxa columns:", len(taxa))
    print("Classified reads:", classified_reads)
    print("Unclassified at rank:", unclassified_reads)
    print(
        "Total represented:",
        classified_reads + unclassified_reads
    )

grand_total = sum(sample_totals.values())

print("\n=== LIB2 ECOLOGY INPUT SUMMARY ===")
print("Sample column:", repr(sample_col))
print("Samples:", len(rows))
print("ASVs:", len(asv_cols))
print("Total LIB2 reads:", grand_total)
print("Expected reads: 186136")
print("Read-total match:", grand_total == 186136)

print("\nCreated ecology input files in:", outdir)
