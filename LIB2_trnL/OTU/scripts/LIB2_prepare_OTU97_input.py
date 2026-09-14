import csv

table_file = "taxonomy/input/LIB2_final_ASV_table.csv"
fasta_file = "taxonomy/input/LIB2_final_ASVs.fasta"
out_fasta = "taxonomy/otu_97/LIB2_ASVs_with_size.fasta"

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

written = 0

with open(fasta_file, encoding="utf-8") as inp, \
     open(out_fasta, "w", encoding="utf-8") as out:

    for line in inp:

        if line.startswith(">"):
            asv = line[1:].strip().split()[0]

            if asv not in totals:
                raise ValueError(
                    f"{asv} missing from abundance table"
                )

            out.write(f">{asv};size={totals[asv]};\n")
            written += 1

        else:
            out.write(line)

total_reads = sum(totals.values())

print("=== LIB2 OTU97 INPUT ===")
print("ASVs in table:", len(asv_cols))
print("FASTA sequences written:", written)
print("Total reads:", total_reads)
print("Expected ASVs: 341")
print("Expected reads: 186136")
print("ASV count match:", len(asv_cols) == 341)
print("FASTA/table match:", written == len(asv_cols))
print("Read-total match:", total_reads == 186136)
print("Output:", out_fasta)
