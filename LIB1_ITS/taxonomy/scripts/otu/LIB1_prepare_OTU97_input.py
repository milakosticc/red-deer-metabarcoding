import csv
import re

table_file = "taxonomy/final/LIB1_target_plant_ASV_table.csv"
fasta_file = "taxonomy/final/LIB1_target_plant_ASVs.fasta"
out_fasta = "taxonomy/otu_97/LIB1_ASVs_with_size.fasta"

# Read abundance table
with open(table_file, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [c for c in headers if c.startswith("ASV")]

totals = {}

for asv in asv_cols:
    total = 0

    for row in rows:
        try:
            total += int(float(row.get(asv, 0) or 0))
        except ValueError:
            pass

    totals[asv] = total

# Read FASTA and write abundance-tagged FASTA
written = 0
current = None

with open(fasta_file, encoding="utf-8") as inp, \
     open(out_fasta, "w", encoding="utf-8") as out:

    for line in inp:

        if line.startswith(">"):
            asv = line[1:].strip().split()[0]

            if asv not in totals:
                raise ValueError(f"{asv} missing from abundance table")

            current = asv
            out.write(f">{asv};size={totals[asv]};\n")
            written += 1

        else:
            out.write(line)

print("=== LIB1 OTU INPUT ===")
print("ASVs:", len(asv_cols))
print("FASTA sequences written:", written)
print("Total reads:", sum(totals.values()))
print("Expected reads: 360007")
print("Read-total match:", sum(totals.values()) == 360007)
print("Output:", out_fasta)
