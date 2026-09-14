import csv
import re

uc_file = "taxonomy/otu_97/LIB1_OTU97_clusters.uc"
asv_table = "taxonomy/final/LIB1_target_plant_ASV_table.csv"
centroids_in = "taxonomy/otu_97/LIB1_OTU97_centroids_raw.fasta"

mapping_out = "taxonomy/otu_97/LIB1_ASV_to_OTU97.tsv"
table_out = "taxonomy/otu_97/LIB1_OTU97_abundance_counts.tsv"
centroids_out = "taxonomy/otu_97/LIB1_OTU97_centroids.fasta"

def clean_id(x):
    return re.sub(r";size=\d+;?$", "", x)

# --------------------------------------------------
# Read UC clustering and make ASV -> OTU mapping
# --------------------------------------------------
cluster_to_otu = {}
asv_to_otu = {}
centroid_to_otu = {}

with open(uc_file, encoding="utf-8") as f:
    for line in f:
        if not line.strip():
            continue

        p = line.rstrip("\n").split("\t")
        record_type = p[0]

        if record_type not in {"S", "H"}:
            continue

        cluster = int(p[1])

        if cluster not in cluster_to_otu:
            cluster_to_otu[cluster] = f"OTU{cluster + 1}"

        otu = cluster_to_otu[cluster]

        query = clean_id(p[8])

        asv_to_otu[query] = otu

        if record_type == "S":
            centroid_to_otu[query] = otu

# --------------------------------------------------
# Write mapping
# --------------------------------------------------
with open(mapping_out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["ASV_ID", "OTU_ID"])

    for asv in sorted(
        asv_to_otu,
        key=lambda x: int(x.replace("ASV", ""))
    ):
        w.writerow([asv, asv_to_otu[asv]])

# --------------------------------------------------
# Read original final ASV count table
# --------------------------------------------------
with open(asv_table, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    headers = reader.fieldnames

asv_cols = [c for c in headers if c in asv_to_otu]
metadata_cols = [c for c in headers if c not in asv_to_otu]

if not metadata_cols:
    raise ValueError("Sample column not detected.")

sample_col = metadata_cols[0]

otus = sorted(
    set(asv_to_otu.values()),
    key=lambda x: int(x.replace("OTU", ""))
)

# --------------------------------------------------
# Build sample x OTU count matrix
# --------------------------------------------------
grand_total = 0

with open(table_out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["Sample"] + otus)

    for row in rows:
        counts = {otu: 0 for otu in otus}

        for asv in asv_cols:
            try:
                n = int(float(row.get(asv, 0) or 0))
            except ValueError:
                n = 0

            counts[asv_to_otu[asv]] += n

        sample_total = sum(counts.values())
        grand_total += sample_total

        w.writerow(
            [row[sample_col]] +
            [counts[otu] for otu in otus]
        )

# --------------------------------------------------
# Rename representative sequences OTU1, OTU2...
# --------------------------------------------------
written = 0
keep = False

with open(centroids_in, encoding="utf-8") as inp, \
     open(centroids_out, "w", encoding="utf-8") as out:

    for line in inp:

        if line.startswith(">"):
            raw = line[1:].strip().split()[0]
            asv = clean_id(raw)

            if asv not in centroid_to_otu:
                raise ValueError(
                    f"Centroid {asv} missing from UC file"
                )

            otu = centroid_to_otu[asv]

            out.write(f">{otu}|centroid={asv}\n")
            written += 1

        else:
            out.write(line)

print("=== LIB1 OTU97 FINAL TABLE VALIDATION ===")
print("Input ASVs mapped:", len(asv_to_otu))
print("Expected ASVs: 466")
print("OTUs:", len(otus))
print("Representative sequences:", written)
print("Samples:", len(rows))
print("Total reads:", grand_total)
print("Expected reads: 360007")
print("Read-total match:", grand_total == 360007)
print("All ASVs mapped:", len(asv_to_otu) == 466)
print("OTU/centroid match:", len(otus) == written)

print("\nCreated:")
print(mapping_out)
print(table_out)
print(centroids_out)
