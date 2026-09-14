import csv
import os
import shutil

counts_in = "taxonomy/otu_97/LIB3_OTU97_abundance_counts.tsv"
taxonomy_in = "taxonomy/otu_97/LIB3_OTU97_taxonomy.tsv"
metadata_in = "taxonomy/ecology_input/LIB3_metadata_template.tsv"

outdir = "taxonomy/ecology_input/OTU97"
os.makedirs(outdir, exist_ok=True)

counts_out = f"{outdir}/LIB3_OTU97_feature_abundance_counts.tsv"
taxonomy_out = f"{outdir}/LIB3_OTU97_taxonomy.tsv"
metadata_out = f"{outdir}/LIB3_OTU97_metadata_template.tsv"

shutil.copyfile(counts_in, counts_out)
shutil.copyfile(taxonomy_in, taxonomy_out)
shutil.copyfile(metadata_in, metadata_out)

with open(counts_out, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    rows = list(reader)
    fields = reader.fieldnames

otu_ids = fields[1:]

total_reads = 0
for row in rows:
    for otu in otu_ids:
        total_reads += int(float(row[otu] or 0))

with open(taxonomy_out, newline="", encoding="utf-8") as f:
    reader = csv.DictReader(f, delimiter="\t")
    tax_rows = list(reader)

taxonomy_ids = {row["OTU_ID"] for row in tax_rows}

print("=== LIB3 OTU97 ECOLOGY INPUT ===")
print("Samples:", len(rows))
print("OTUs:", len(otu_ids))
print("Taxonomy rows:", len(tax_rows))
print("Total reads:", total_reads)
print("Expected reads: 307096")
print("Read-total match:", total_reads == 307096)
print("OTU IDs match taxonomy:", set(otu_ids) == taxonomy_ids)

print("\nFiles:")
print(counts_out)
print(taxonomy_out)
print(metadata_out)
