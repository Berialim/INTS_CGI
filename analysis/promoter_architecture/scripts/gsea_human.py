# -*- coding: utf-8 -*-
"""
This script performs Gene Set Enrichment Analysis (GSEA) on a differential
expression results file. It includes a preprocessing step to map RefSeq IDs 
(NM_...) to Gene Symbols (name2) using a reference file.

Dependencies:
pip install pandas gseapy matplotlib
"""

import pandas as pd
import gseapy as gp
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from gseapy import dotplot

# --- Configuration ---
INPUT_FILE = "results.txt"
REFSEQ_FILE = "RefSeq_all.txt"
GENE_SET_LIBRARIES = {
    'GO_BP': 'GO_Biological_Process_2023',
    #'GO_MF': 'GO_Molecular_Function_2023',
    'KEGG': 'KEGG_2021_Human',
    'HALLMARK': 'h.all.v2025.1.Hs.symbols.gmt' 
}

def map_refseq_to_symbol(df, ref_path):
    """
    Checks if the index contains RefSeq IDs (NM_) and maps them to Gene Symbols (name2).
    """
    # Clean the index (remove whitespace, quotes)
    df.index = df.index.astype(str).str.strip('"').str.strip("'")
    first_id = df.index[0]
    
    if not (first_id.startswith('NM_') or first_id.startswith('NR_')):
        print(f"First ID '{first_id}' does not look like RefSeq. Skipping mapping...")
        return df

    print(f"RefSeq IDs detected. Mapping using {ref_path}...")
    
    try:
        # Load reference file. Note: some RefSeq files have a leading # on 'bin' column
        ref_df = pd.read_csv(os.path.expanduser(ref_path), sep='\t')
        
        # Ensure column names are clean (remove leading # if present)
        ref_df.columns = [c.lstrip('#') for c in ref_df.columns]
        
        # Create mapping dictionaries
        full_map = dict(zip(ref_df['name'].astype(str), ref_df['name2'].astype(str)))
        
        # Also create a map for stripped versions (NM_123.1 -> NM_123)
        ref_df['name_no_ver'] = ref_df['name'].astype(str).str.split('.').str[0]
        stripped_map = dict(zip(ref_df['name_no_ver'], ref_df['name2'].astype(str)))
        
    except Exception as e:
        print(f"Error loading reference file: {e}")
        return df

    # Prepare IDs for matching
    original_ids = df.index.tolist()
    stripped_ids = [str(x).split('.')[0] for x in original_ids]
    
    # Attempt Mapping
    symbols = [full_map.get(oid) for oid in original_ids]
    
    # Fill in gaps using stripped IDs
    for i, sym in enumerate(symbols):
        if pd.isna(sym):
            symbols[i] = stripped_map.get(stripped_ids[i])
            
    df['Gene_Symbol'] = symbols
    
    # IMPORTANT: Drop rows that failed to map BEFORE deduplication
    unmapped_count = df['Gene_Symbol'].isna().sum()
    df = df.dropna(subset=['Gene_Symbol'])
    
    if unmapped_count > 0:
        print(f"Warning: {unmapped_count} IDs could not be mapped and were removed.")

    if df.empty:
        print("DEBUG: Final DataFrame is empty after mapping.")
        print(f"Sample input ID from your file: '{original_ids[0]}'")
        print(f"Sample Ref ID from your reference: '{list(full_map.keys())[0]}'")
        return df

    # Handle duplicates: Multiple NM_ mapping to one Gene Symbol.
    # We sort by absolute log2FC and keep the entry with the highest magnitude for each gene.
    df['abs_fc'] = df['log2FoldChange'].abs()
    df = df.sort_values('abs_fc', ascending=False).drop_duplicates('Gene_Symbol')
    
    df = df.set_index('Gene_Symbol')
    print(f"Mapping complete. Final gene count: {len(df)}")
    return df

def run_gsea():
    print(f"Loading data from '{INPUT_FILE}'...")
    if not os.path.exists(INPUT_FILE):
        print(f"Error: The file '{INPUT_FILE}' was not found.")
        sys.exit(1)

    try:
        df = pd.read_csv(INPUT_FILE, index_col=0)
        df.index.name = 'Gene_Name'
    except Exception as e:
        print(f"Error reading the file: {e}")
        sys.exit(1)

    # 1. Map IDs
    df = map_refseq_to_symbol(df, REFSEQ_FILE)

    # 2. Preprocessing & Ranking
    # Crucial: Filter for necessary columns
    df_cleaned = df.dropna(subset=['log2FoldChange', 'padj']).copy()
    
    if df_cleaned.empty:
        print("Error: No valid data found after cleaning. Check log2FoldChange/padj column names.")
        sys.exit(1)
    
    df_cleaned.dropna(inplace=True)
    # Prevent log10(0)
    df_cleaned['padj'] = df_cleaned['padj'].replace(0, 2e-302)
    
    # Calculate ranking score: sign(log2FC) * -log10(padj)
    df_cleaned['ranking_score'] = df_cleaned['log2FoldChange']  # * -np.log10(df_cleaned['padj'])
    
    ranked_genes = df_cleaned.sort_values(by='ranking_score', ascending=False)
    
    # 3. GSEA Loop
    for name, library in GENE_SET_LIBRARIES.items():
        output_dir = f'gsea_results_{name}'
        print(f"\n--- Running GSEA for '{name}' ---")
        
        try:
            prerank_results = gp.prerank(
                rnk=ranked_genes['ranking_score'],
                gene_sets=library,
                outdir=output_dir,
                permutation_num=1000,
                min_size=5,
                max_size=1000,
                threads=18,
                no_plot=False
            )

            res = prerank_results.res2d
            if res.empty:
                print(f"No results found for {name}.")
                continue
                
            terms = res['Term'].tolist()
            
            # Specific pathway plot
            pi3k_terms = [t for t in terms if 'PI3K' in t or 'AKT' in t.upper()]
            if pi3k_terms:
                prerank_results.plot(terms=pi3k_terms[0:1], ofname=os.path.join(output_dir, "PI3K_Akt_GSEA_plot.pdf"))
            
            # Dotplot
            dotplot(res, 
                    ofname=os.path.join(output_dir, "GSEA_dotplot.pdf"), 
                    column="FDR q-val",
                    title=f'{name} Enrichment',
                    cmap=plt.cm.RdBu_r,
                    cutoff=0.25, 
                    figsize=(6, 8))

            # Top 5 plot
            if len(terms) >= 5:
                prerank_results.plot(terms[:5], show_ranking=True, ofname=os.path.join(output_dir, "GSEA_top5plot.pdf"))

            print(f"Done. Results saved in {output_dir}")

        except Exception as e:
            print(f"Error running GSEA for {name}: {e}")

if __name__ == "__main__":
    run_gsea()
