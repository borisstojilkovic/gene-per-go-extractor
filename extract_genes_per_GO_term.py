# -*- coding: utf-8 -*-
"""
Extract genes per GO term from RNA-seq result files, annotate them, and export:
1) a per-GO table with merged RNA-seq + annotations
2) a grouped table where each column is a GO term and rows are GeneIDs

Inputs:
- "input/" folder with RNA-seq outputs (DESeq2 .tab or .xlsx/.xls)
- "Go_termnIDs_and_file_names.xlsx" listing GO terms in column "names"
- GO legend files (TSV) and annotation Excel files in "annotations/" and project root

Outputs:
- "output/<basename>/" per GO tables and a grouped table

Usage:
- Run the script, choose species in the prompt (S, SW, or A).
"""

import os
import csv
import math
import pandas as pd

# pd.set_option("display.max_rows", None, "display.max_columns", None)

print(
    "#######\n"
    "Welcome to the Gene-per-GO extractor\n\n"
    "• Open: Go_termnIDs_and_file_names.xlsx\n"
    "  - Put the GO terms you want to process in the first column (e.g., GO:0009535)\n"
    "• Place all RNA-seq result files in the 'input' folder.\n"
    "• Outputs will be written to 'output/<file_basename>/'.\n"
    "########"
)

# Ensure the output folder exists
if not os.path.exists("output"):
    os.mkdir("output")

# Species selection
species = input(
    "Select species:\n"
    "  - Type 'S' for tomato\n"
    "  - Type 'SW' for tomato with RKN (Mi–Tomato)\n"
    "  - Type 'A' for Arabidopsis\n"
    "Your choice: "
).strip().upper()

# Load species-specific annotation and GO legend
if species == "S":
    annotation_file = pd.read_excel("annotations/annotation_tom.xlsx")
    go_legend_path = "GO term accession_mart_exportSL_3.0_with descr.tab"
elif species == "SW":
    annotation_file = pd.read_excel("annotations/annotation_tom_with_RKN_Mi-Tomato.xlsx")
    go_legend_path = "GO term accession_mart_exportSL_3.0_with descr.tab"
elif species == "A":
    annotation_file = pd.read_excel("annotations/annotation_arab.xlsx")
    go_legend_path = "GO term accession_mart_exportArabidopsis.tab"
else:
    print("Unrecognized species option. Please run again and choose S, SW, or A.")
    raise SystemExit(1)

try:
    Go_legend = pd.read_csv(go_legend_path, sep="\t", low_memory=False)
except FileNotFoundError:
    print(f"GO legend file not found: {go_legend_path}")
    raise

# Load GO terms to process
try:
    Go_termnIDs = pd.read_excel("Go_termnIDs_and_file_names.xlsx")
except FileNotFoundError:
    print("File 'Go_termnIDs_and_file_names.xlsx' not found. Please add it and try again.")
    raise

if "names" not in Go_termnIDs.columns:
    print("Expected a column named 'names' in Go_termnIDs_and_file_names.xlsx.")
    raise SystemExit(1)

Go_termnIDs_list = Go_termnIDs["names"].tolist()

# Process each RNA-seq file in the input folder
for file in os.listdir("input"):
    # Skip non-files (e.g., directories)
    in_path = os.path.join("input", file)
    if not os.path.isfile(in_path):
        continue

    # Create an output subfolder per input file (by basename)
    base, ext = os.path.splitext(file)
    sub_out = os.path.join("output", base)
    if not os.path.exists(sub_out):
        os.mkdir(sub_out)

    # Load RNA-seq table
    if ext.lower() == ".tab":
        RNA_seqfile = pd.read_csv(in_path, sep="\t", low_memory=False)
    elif ext.lower() in (".xlsx", ".xls"):
        RNA_seqfile = pd.read_excel(in_path)
    else:
        print(f"Skipping unsupported file type: {file}")
        continue

    print(f"Processing: {file}")

    # Ensure GeneID is present
    if "GeneID" not in RNA_seqfile.columns:
        print(f"  Warning: 'GeneID' column not found in {file}. Skipping this file.")
        continue

    # Prepare locus column for merging
    if species == "A":
        RNA_seqfile["locus"] = RNA_seqfile["GeneID"]
    else:
        # For tomato-based IDs, split at the first dot to get the locus
        split_cols = RNA_seqfile["GeneID"].astype(str).str.split(".", n=1, expand=True)
        RNA_seqfile["locus"] = split_cols[0]

    # Prepare holder for the grouped table (columns: GO terms; rows: GeneIDs)
    grouped_df = pd.DataFrame({"_pad": list(range(10000))})  # padding to equalize column lengths

    # Loop over GO terms
    for go in Go_termnIDs_list:
        if isinstance(go, float) and math.isnan(go):
            continue

        go = str(go).strip()
        go_safe = go.replace(":", "_")

        # Subset GO legend for this accession
        df_go = Go_legend.loc[Go_legend["accession"].isin([go]), ["Gene stable ID", "GO term name", "GO term definition"]]

        if df_go.empty:
            print(f"  Note: No entries found for GO term {go}. Skipping.")
            continue

        # Get GO term name (safe fallback)
        term_name = df_go["GO term name"].iloc[0] if "GO term name" in df_go.columns else go
        term_name_safe = str(term_name).replace("/", "_")

        # Save raw GO legend subset for this term
        go_legend_out = os.path.join(sub_out, f"{go_safe} {term_name_safe}.tab")
        df_go.to_csv(go_legend_out, sep="\t", index=False)

        # Reload (optional, keeps your original flow)
        df = pd.read_csv(go_legend_out, sep="\t", low_memory=False)

        # Create locus for merging
        if species == "A":
            df["locus"] = df["Gene stable ID"]
        else:
            df[["locus", "_rest"]] = df["Gene stable ID"].astype(str).str.split(".", n=1, expand=True)

        # Merge with RNA-seq and annotations
        df = pd.merge(df, RNA_seqfile, on="locus", how="inner")
        df = pd.merge(df, annotation_file, on="locus", how="left")

        # Clean up temporary columns
        cols_to_drop = [c for c in ["_rest", "locus"] if c in df.columns]
        if "Gene stable ID" in df.columns:
            # Remove duplicate rows by the stable ID
            df = df.drop_duplicates(subset=["Gene stable ID"])
            cols_to_drop.append("Gene stable ID")

        if cols_to_drop:
            df = df.drop(columns=cols_to_drop, errors="ignore")

        # Reorder columns with GeneID first (if present)
        if "GeneID" in df.columns:
            other_cols = [c for c in df.columns if c != "GeneID"]
            df = df[["GeneID"] + other_cols]

        # Save per-GO merged table
        per_go_out = os.path.join(sub_out, f"{base}{go_safe} {term_name_safe}.tab")
        df.to_csv(per_go_out, sep="\t", index=False, decimal=",")

        # Build grouped table column (GeneID list, padded to 10000)
        if "GeneID" in df.columns:
            geneids = df["GeneID"].astype(str).tolist()
        else:
            # Fallback: if GeneID missing, try to use any gene-like column
            geneids = []

        if len(geneids) < 10000:
            geneids.extend([""] * (10000 - len(geneids)))
        else:
            geneids = geneids[:10000]

        grouped_df[term_name_safe] = geneids

    # Finalize grouped table
    if "_pad" in grouped_df.columns:
        grouped_df = grouped_df.drop(columns=["_pad"])

    grouped_out = os.path.join(sub_out, f"{base}1_grouped.tab")
    grouped_df.to_csv(grouped_out, sep="\t", index=False, decimal=",")

print("Finished!")
