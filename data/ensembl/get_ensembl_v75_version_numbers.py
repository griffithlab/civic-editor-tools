#!/usr/bin/env python3

import pymysql
import csv
import os
import sys

#set data path location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(SCRIPT_DIR, "build37", "ensembl75_transcripts.tsv")

if os.path.exists(output_path):
    print(f"{output_path} already exists")
    sys.exit(0)

# Ensembl public MySQL — no password needed
conn = pymysql.connect(
    host="ensembldb.ensembl.org",
    port=5306,
    user="anonymous",
    password="",
    database="homo_sapiens_core_75_37",  # Ensembl 75 / GRCh37
)

query = """
SELECT
    t.stable_id        AS transcript_id,
    t.version          AS transcript_version,
    g.stable_id        AS gene_id,
    g.version          AS gene_version,
    xref.display_label AS gene_name,
    t.biotype
FROM transcript t
JOIN gene g ON t.gene_id = g.gene_id
LEFT JOIN xref ON g.display_xref_id = xref.xref_id
"""

with conn.cursor() as cur:
    cur.execute(query)
    rows = cur.fetchall()
    columns = [d[0] for d in cur.description]

conn.close()

os.makedirs(os.path.dirname(output_path), exist_ok=True)

with open(output_path, "w", newline="") as f:
    writer = csv.writer(f, delimiter="\t")
    writer.writerow(columns)
    writer.writerows(rows)

print(f"Done. {len(rows):,} transcripts written.")
