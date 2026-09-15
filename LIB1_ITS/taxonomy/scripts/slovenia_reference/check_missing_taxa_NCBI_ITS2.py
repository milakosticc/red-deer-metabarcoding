#!/usr/bin/env python3

import argparse
import csv
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


parser = argparse.ArgumentParser()

parser.add_argument(
    "--email",
    required=True,
    help="Your email address required/recommended for NCBI E-utilities."
)

parser.add_argument(
    "--api-key",
    default=None,
    help="Optional NCBI API key."
)

args = parser.parse_args()


# --------------------------------------------------
# Files
# --------------------------------------------------

input_file = Path(
    "taxonomy/reference/Slovenia_ITS2/missing_species.txt"
)

output_file = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "Slovenia_missing_taxa_NCBI_ITS2_hits.tsv"
)

log_file = Path(
    "taxonomy/reference/Slovenia_ITS2/NCBI/"
    "NCBI_ITS2_search.log"
)


# --------------------------------------------------
# NCBI settings
# --------------------------------------------------

base_url = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
)

# Stay safely below NCBI limits.
delay = 0.38 if not args.api_key else 0.12


# --------------------------------------------------
# Read taxa
# --------------------------------------------------

with open(input_file) as f:

    species_list = [
        line.strip()
        for line in f
        if line.strip()
    ]


print(f"Taxa to search: {len(species_list)}")


# --------------------------------------------------
# Search function
# --------------------------------------------------

def search_ncbi(species):

    # Example:
    #
    # "Achillea clavenae"[Organism]
    # AND
    # ("internal transcribed spacer 2"[Title] OR ITS2[Title])

    query = (
        f'"{species}"[Organism] AND '
        f'("internal transcribed spacer 2"[Title] OR ITS2[Title])'
    )

    params = {
        "db": "nuccore",
        "term": query,
        "retmode": "xml",
        "retmax": "0",
        "tool": "red_deer_metabarcoding",
        "email": args.email,
    }

    if args.api_key:
        params["api_key"] = args.api_key

    url = base_url + "?" + urllib.parse.urlencode(params)

    for attempt in range(1, 4):

        try:

            with urllib.request.urlopen(
                url,
                timeout=30
            ) as response:

                xml_data = response.read()

            root = ET.fromstring(xml_data)

            count = int(
                root.findtext("Count", default="0")
            )

            return count, query, "OK"

        except Exception as e:

            if attempt == 3:
                return 0, query, f"ERROR: {e}"

            time.sleep(2 * attempt)


# --------------------------------------------------
# Search all missing taxa
# --------------------------------------------------

results = []

with open(log_file, "w") as log:

    for i, species in enumerate(species_list, start=1):

        count, query, status = search_ncbi(species)

        results.append({
            "Species": species,
            "NCBI_ITS2_hits": count,
            "Status": status
        })

        message = (
            f"[{i}/{len(species_list)}] "
            f"{species}: {count} ITS2 hits"
        )

        print(message, flush=True)
        log.write(message + "\n")
        log.flush()

        time.sleep(delay)


# --------------------------------------------------
# Save table
# --------------------------------------------------

with open(output_file, "w", newline="") as f:

    writer = csv.DictWriter(
        f,
        fieldnames=[
            "Species",
            "NCBI_ITS2_hits",
            "Status"
        ],
        delimiter="\t",
        lineterminator="\n"
    )

    writer.writeheader()
    writer.writerows(results)


taxa_with_hits = sum(
    1
    for r in results
    if r["NCBI_ITS2_hits"] > 0
)

taxa_without_hits = sum(
    1
    for r in results
    if r["NCBI_ITS2_hits"] == 0
)

total_hits = sum(
    r["NCBI_ITS2_hits"]
    for r in results
)


print()
print("=== NCBI ITS2 SEARCH COMPLETE ===")
print(f"Taxa searched:             {len(results)}")
print(f"Taxa with ITS2 hits:       {taxa_with_hits}")
print(f"Taxa without ITS2 hits:    {taxa_without_hits}")
print(f"Total NCBI ITS2 records:   {total_hits}")
print()
print(f"Results: {output_file}")
print(f"Log:     {log_file}")
