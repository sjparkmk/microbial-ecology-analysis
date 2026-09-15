\
\
\
\
\
\
\
\

import csv

IN_CSV = "/input/youre/data/asv_taxonomy_with_bootstrap_ranked.csv"
OUT_FASTA = "/output/your/data/all_663_asvs.fasta"
OUT_META = "/output/your/data/asv_metadata.tsv"

import os
os.makedirs("/output/your/data", exist_ok=True)

rows = []
with open(IN_CSV, newline="") as f:
    r = csv.DictReader(f)
    for i, row in enumerate(r, start=1):
        asv_id = f"ASV_{i}"
        rows.append((asv_id, row))

with open(OUT_FASTA, "w") as fasta, open(OUT_META, "w", newline="") as meta_f:
    fasta_lines = []
    meta_writer = csv.writer(meta_f, delimiter="\t")
    meta_writer.writerow([
        "ASV_ID", "rank", "sequence", "total_reads", "pct_of_total",
        "Kingdom", "Kingdom_boot", "Superphylum", "Superphylum_boot",
        "Phylum", "Phylum_boot", "Class", "Class_boot", "Order", "Order_boot",
        "Family", "Family_boot", "Genus", "Genus_boot", "Species", "Species_boot",
    ])
    for i, (asv_id, row) in enumerate(rows, start=1):
        seq = row["ASV"]
        fasta_lines.append(f">{asv_id}\n{seq}\n")
        meta_writer.writerow([
            asv_id, i, seq, row["total_reads"], row["pct_of_total"],
            row["Kingdom"], row["Kingdom_boot"], row["Superphylum"], row["Superphylum_boot"],
            row["Phylum"], row["Phylum_boot"], row["Class"], row["Class_boot"],
            row["Order"], row["Order_boot"], row["Family"], row["Family_boot"],
            row["Genus"], row["Genus_boot"], row["Species"], row["Species_boot"],
        ])
    fasta.writelines(fasta_lines)

print(f"Wrote {len(rows)} ASV sequences to {OUT_FASTA}")
print(f"Wrote metadata to {OUT_META}")
