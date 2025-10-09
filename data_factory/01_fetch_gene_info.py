import pandas as pd
import os

# --- Configuration ---
# Define the output directory and filename
OUTPUT_DIR = "data_factory/output"
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "human_gene_info.csv")

# Ensembl BioMart URL and query parameters
BIOMART_URL = "http://www.ensembl.org/biomart/martservice"
DATASET = "hsapiens_gene_ensembl"
ATTRIBUTES = [
    "ensembl_gene_id",
    "entrezgene_id",
    "hgnc_symbol"
]

print("--- Starting Step 1: Fetch Human Gene Info ---")

# --- Ensure output directory exists ---
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)
    print(f"Created output directory: {OUTPUT_DIR}")

# --- Construct the BioMart Query XML ---
query_xml = f"""
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "CSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
    <Dataset name = "{DATASET}" interface = "default" >
        <Attribute name = "{ATTRIBUTES[0]}" />
        <Attribute name = "{ATTRIBUTES[1]}" />
        <Attribute name = "{ATTRIBUTES[2]}" />
    </Dataset>
</Query>
"""

try:
    # --- Fetch data from Ensembl BioMart ---
    print(f"Connecting to Ensembl BioMart to fetch dataset: {DATASET}...")
    
    # Use pandas to read the data directly from the URL
    # Pandas can handle the web request and parse the CSV response
    df = pd.read_csv(BIOMART_URL, params={'query': query_xml})
    
    # Assign column names since header=0
    df.columns = ATTRIBUTES
    
    print(f"Successfully fetched data for {len(df)} genes.")
    
    # --- Data Cleaning ---
    # Drop rows where essential identifiers are missing
    df.dropna(subset=['ensembl_gene_id', 'hgnc_symbol'], inplace=True)
    
    # Convert Entrez ID to integer, handling potential missing values
    # Coerce errors will turn non-numeric values into NaT (Not a Time) which we then fill
    df['entrezgene_id'] = pd.to_numeric(df['entrezgene_id'], errors='coerce').fillna(0).astype(int)
    
    print(f"Cleaned data, resulting in {len(df)} genes with complete core info.")

    # --- Save to CSV ---
    df.to_csv(OUTPUT_FILE, index=False)
    
    print(f"Successfully saved gene information to: {OUTPUT_FILE}")
    print("--- Step 1: COMPLETED ---")

except Exception as e:
    print(f"An error occurred: {e}")
    print("--- Step 1: FAILED ---")
