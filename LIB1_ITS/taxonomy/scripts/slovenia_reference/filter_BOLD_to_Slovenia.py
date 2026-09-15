#!/usr/bin/env python3

from pathlib import Path

base = Path("taxonomy/reference")

whitelist_file = base / "Slovenia_ITS2/Slovenia_plant_species.txt"
bold_file = base / "BOLD_plant_ITS/BOLD_Plantae_ITS2_SINTAX.fasta"

output_fasta = base / "Slovenia_ITS2/Slovenia_BOLD_Plantae_ITS2_SINTAX.fasta"
matched_file = base / "Slovenia_ITS2/matched_species.txt"
missing_file = base / "Slovenia_ITS2/missing_species.txt"


# -----------------------------
# Read Slovenian species list
# -----------------------------

with open(whitelist_file) as f:
    slovenia_species = {
        line.strip()
        for line in f
        if line.strip()
    }


# -----------------------------
# Read BOLD FASTA
# -----------------------------

records_kept = 0
bold_species_found = set()

keep_record = False

with open(bold_file) as infile, open(output_fasta, "w") as outfile:

    for line in infile:

        if line.startswith(">"):

            keep_record = False
            species = None

            # Example:
            # s:Codiaeum_variegatum;
            if "s:" in line:
                species_part = line.split("s:", 1)[1]
                species = species_part.split(";", 1)[0]
                species = species.replace("_", " ").strip()

            if species and species in slovenia_species:
                keep_record = True
                records_kept += 1
                bold_species_found.add(species)

        if keep_record:
            outfile.write(line)


# -----------------------------
# Match statistics
# -----------------------------

matched_species = sorted(slovenia_species & bold_species_found)
missing_species = sorted(slovenia_species - bold_species_found)


with open(matched_file, "w") as f:
    for species in matched_species:
        f.write(species + "\n")


with open(missing_file, "w") as f:
    for species in missing_species:
        f.write(species + "\n")


print("=== Slovenia ITS2 database filtering ===")
print(f"Species in Slovenia whitelist: {len(slovenia_species)}")
print(f"Species matched in BOLD ITS2:   {len(matched_species)}")
print(f"Species missing from BOLD ITS2: {len(missing_species)}")
print(f"BOLD reference sequences kept:  {records_kept}")
print()
print(f"Output FASTA: {output_fasta}")
print(f"Matched list: {matched_file}")
print(f"Missing list: {missing_file}")
