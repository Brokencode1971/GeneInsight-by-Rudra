import pandas as pd
import os
import requests

# --- HELPER FUNCTION TO PARSE GO DICTIONARY ---
def parse_obo(file_path):
    go_terms = []
    with open(file_path, 'r', encoding='utf-8') as f:
        term_data = {}
        for line in f:
            line = line.strip()
            if line == '[Term]':
                if term_data and 'id' in term_data:
                    go_terms.append({'GO_ID': term_data['id'], 'GO_Term': term_data.get('name', '')})
                term_data = {}
            elif ':' in line:
                key, value = line.split(':', 1)
                term_data[key.strip()] = value.strip()
    if term_data and 'id' in term_data:
        go_terms.append({'GO_ID': term_data['id'], 'GO_Term': term_data.get('name', '')})
    return pd.DataFrame(go_terms)

# --- SETUP ---
RAW_DATA_DIR = 'raw_data'
PROCESSED_DATA_DIR = 'processed_data'
if not os.path.exists(PROCESSED_DATA_DIR):
    os.makedirs(PROCESSED_DATA_DIR)

def raw_path(filename):
    return os.path.join(RAW_DATA_DIR, filename)

def processed_path(filename):
    return os.path.join(PROCESSED_DATA_DIR, filename)

# --- 0. DOWNLOAD AND PROCESS GO DICTIONARY ---
print("Processing GO dictionary (go-basic.obo)...")
obo_file_path = raw_path('go-basic.obo')
# (Assuming this file was downloaded by a previous script version)
if os.path.exists(obo_file_path):
    go_map_df = parse_obo(obo_file_path)
    go_map_df.to_csv(processed_path('go_terms.tsv'), sep='\t', index=False)
    print("-> Clean 'go_terms.tsv' created.")
else:
    print("-> WARNING: go-basic.obo not found. Skipping GO term processing.")


# --- 1. PROCESS ENSEMBL TO GO MAPPING ---
print("Processing Ensembl to GO mapping file...")
ensembl_go_path = raw_path('mart_export.txt.gz') # This is the new file from BioMart
if os.path.exists(ensembl_go_path):
    ensembl_go_df = pd.read_csv(ensembl_go_path, compression='gzip', sep='\t')
    # Rename columns for clarity
    ensembl_go_df.rename(columns={
        'Gene stable ID': 'ensembl_gene_id',
        'GO term accession': 'go_id'
    }, inplace=True)
    # Drop any rows where GO ID is missing
    ensembl_go_df.dropna(subset=['go_id'], inplace=True)
    ensembl_go_df.to_csv(processed_path('ensembl_to_go.tsv'), sep='\t', index=False)
    print("-> Clean 'ensembl_to_go.tsv' saved.")
else:
     print(f"-> WARNING: {os.path.basename(ensembl_go_path)} not found. Skipping.")


# --- 2. PROCESS OTHER FILES (GENE INFO, PPIs, etc.) ---
# This part of the script remains the same as before, processing the other files.
print("Processing remaining data files...")

# Process Gene Info
try:
    gene_info_df = pd.read_csv(raw_path('mart_export.txt.gz'), compression='gzip', sep='\t')
    gene_info_df.rename(columns={
        'Gene stable ID': 'ensembl_gene_id', 'Gene name': 'gene_symbol'
    }, inplace=True)
    # Select only the columns we need for mapping
    gene_map_df = gene_info_df[['ensembl_gene_id', 'gene_symbol']].dropna().drop_duplicates()
    gene_map_df.to_csv(processed_path('gene_map.tsv'), sep='\t', index=False)
    print("-> Clean 'gene_map.tsv' saved.")
except FileNotFoundError:
    print("-> WARNING: Could not find original mart_export.txt.gz for gene mapping. Skipping.")


# Process BioGRID
try:
    biogrid_df = pd.read_csv(raw_path('BIOGRID-ORGANISM-Homo_sapiens-5.0.250.tab3.txt'), sep='\t', low_memory=False)
    filtered_biogrid_df = biogrid_df[biogrid_df['Experimental System Type'] == 'physical'].copy()
    filtered_biogrid_df.to_csv(processed_path('filtered_biogrid.tsv'), sep='\t', index=False)
    print("-> Clean 'filtered_biogrid.tsv' saved.")
except FileNotFoundError:
    print("-> WARNING: BioGRID file not found. Skipping.")


print("\n--- SCRIPT COMPLETE ---")

