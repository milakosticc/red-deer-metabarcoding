import csv
from collections import defaultdict

base = "taxonomy/reference/trnL_CH_2026/trnL_DB_Zenodo"

fasta_in = f"{base}/DBs_fasta/Obitools_trnl_CH.fasta"
tax_in   = f"{base}/DBs_taxonomy/Obitools_trnl_CH.csv"

fasta_out = "taxonomy/reference/trnL_CH_2026/Obitools_trnl_CH_SINTAX.fasta"

taxonomy = {}

with open(tax_in, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        taxonomy[row["Acc"]] = row

records = []
current = None
seq = []

with open(fasta_in, encoding="utf-8") as f:
    for line in f:
        line = line.strip()

        if line.startswith(">"):
            if current is not None:
                records.append((current, "".join(seq)))

            current = line[1:].split()[0]
            seq = []
        else:
            seq.append(line)

    if current is not None:
        records.append((current, "".join(seq)))

counts = defaultdict(int)
written = 0
missing = 0

def clean(x):
    x = (x or "").strip()
    if not x or x.lower() in {"na", "nan", "none"}:
        return ""
    return x.replace(" ", "_").replace(";", "_").replace(",", "_")

with open(fasta_out, "w", encoding="utf-8") as out:

    for acc, sequence in records:

        if acc not in taxonomy:
            missing += 1
            continue

        counts[acc] += 1

        unique_id = acc
        if counts[acc] > 1:
            unique_id = f"{acc}_{counts[acc]}"

        t = taxonomy[acc]

        ranks = [
            ("k", clean(t["Kingdom"])),
            ("p", clean(t["Phylum"])),
            ("c", clean(t["Class"])),
            ("o", clean(t["Order"])),
            ("f", clean(t["Family"])),
            ("g", clean(t["Genus"])),
            ("s", clean(t["Species"]))
        ]

        tax_string = ",".join(
            f"{prefix}:{value}"
            for prefix, value in ranks
            if value
        )

        out.write(f">{unique_id};tax={tax_string};\n")
        out.write(sequence.upper() + "\n")

        written += 1

print("Input FASTA records:", len(records))
print("Reference records written:", written)
print("Records without taxonomy:", missing)
print("Output:", fasta_out)
