import pandas as pd
import os
from fastapi import FastAPI
from pydantic import BaseModel
from contextlib import asynccontextmanager
import re
# **NEW:** Import the CORSMiddleware
from fastapi.middleware.cors import CORSMiddleware

# --- DATA STORAGE ---
# Global variables to hold our data in memory
gene_info_df = None
go_terms_map_df = None
biogrid_ppi_df = None
ensembl_to_go_df = None

# --- SETUP: LOAD DATA ON STARTUP ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Server starting up...")
    print("Loading data files into memory...")

    global gene_info_df, go_terms_map_df, biogrid_ppi_df, ensembl_to_go_df
    
    processed_dir = 'processed_data'
    
    gene_info_df = pd.read_csv(os.path.join(processed_dir, 'genes.tsv'), sep='\t')
    go_terms_map_df = pd.read_csv(os.path.join(processed_dir, 'go_terms.tsv'), sep='\t') 
    biogrid_ppi_df = pd.read_csv(os.path.join(processed_dir, 'biogrid_ppi.tsv'), sep='\t', low_memory=False)
    ensembl_to_go_df = pd.read_csv(os.path.join(processed_dir, 'ensembl_to_go.tsv'), sep='\t')

    # Data cleaning and uppercase conversion for consistent matching
    if 'gene_symbol' in gene_info_df.columns:
        gene_info_df['gene_symbol'] = gene_info_df['gene_symbol'].str.upper()
    if 'Official Symbol Interactor A' in biogrid_ppi_df.columns:
        biogrid_ppi_df['Official Symbol Interactor A'] = biogrid_ppi_df['Official Symbol Interactor A'].str.upper()
    if 'Official Symbol Interactor B' in biogrid_ppi_df.columns:
        biogrid_ppi_df['Official Symbol Interactor B'] = biogrid_ppi_df['Official Symbol Interactor B'].str.upper()

    print("Data loading complete. Server is ready.")
    yield
    print("Server shutting down...")

# --- INITIALIZE FASTAPI APP ---
app = FastAPI(lifespan=lifespan)

# **NEW:** Add the CORS middleware to the application
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Allows all origins
    allow_credentials=True,
    allow_methods=["*"],  # Allows all methods
    allow_headers=["*"],  # Allows all headers
)

# --- DEFINE INPUT MODEL ---
class GeneLists(BaseModel):
    up_regulated: list[str]
    down_regulated: list[str]

# --- API ENDPOINT ---
@app.post("/compare")
async def compare_gene_lists(lists: GeneLists):
    # 1. Clean and validate the input Ensembl ID lists
    ensembl_ids_a = {s.strip().upper() for s in lists.up_regulated if s.strip()}
    ensembl_ids_b = {s.strip().upper() for s in lists.down_regulated if s.strip()}

    # 2. Get the count of mapped genes
    db_ensembl_ids = set(gene_info_df['ensembl_gene_id'].dropna())
    mapped_a_count = len(ensembl_ids_a.intersection(db_ensembl_ids))
    mapped_b_count = len(ensembl_ids_b.intersection(db_ensembl_ids))

    # 3. Perform GO Term Comparison
    go_a = ensembl_to_go_df[ensembl_to_go_df['ensembl_gene_id'].isin(ensembl_ids_a)]
    go_b = ensembl_to_go_df[ensembl_to_go_df['ensembl_gene_id'].isin(ensembl_ids_b)]
    
    go_set_a = set(go_a['go_id'])
    go_set_b = set(go_b['go_id'])
    
    unique_go_a_ids = list(go_set_a - go_set_b)
    unique_go_b_ids = list(go_set_b - go_set_a)
    shared_go_ids = list(go_set_a & go_set_b)

    def get_terms_from_ids(id_list):
        if not id_list: return []
        terms = go_terms_map_df[go_terms_map_df['GO_ID'].isin(id_list)]
        return terms.rename(columns={'GO_ID': 'id', 'GO_Term': 'term'}).to_dict('records')

    # 4. Perform PPI Network Analysis
    symbols_a = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_a)]['gene_symbol'])
    symbols_b = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_b)]['gene_symbol'])

    internal_a = biogrid_ppi_df[biogrid_ppi_df['Official Symbol Interactor A'].isin(symbols_a) & biogrid_ppi_df['Official Symbol Interactor B'].isin(symbols_a)]
    internal_b = biogrid_ppi_df[biogrid_ppi_df['Official Symbol Interactor A'].isin(symbols_b) & biogrid_ppi_df['Official Symbol Interactor B'].isin(symbols_b)]
    cross_talk = biogrid_ppi_df[
        (biogrid_ppi_df['Official Symbol Interactor A'].isin(symbols_a) & biogrid_ppi_df['Official Symbol Interactor B'].isin(symbols_b)) |
        (biogrid_ppi_df['Official Symbol Interactor A'].isin(symbols_b) & biogrid_ppi_df['Official Symbol Interactor B'].isin(symbols_a))
    ]
    ppi_cols = ['Official Symbol Interactor A', 'Official Symbol Interactor B']
    
    # 5. Final JSON Response
    return {
        "summary": {
            "up_regulated_submitted_count": len(ensembl_ids_a),
            "down_regulated_submitted_count": len(ensembl_ids_b),
            "up_regulated_mapped_count": mapped_a_count,
            "down_regulated_mapped_count": mapped_b_count,
        },
        "go_comparison": {
            "unique_to_up_regulated": get_terms_from_ids(unique_go_a_ids),
            "unique_to_down_regulated": get_terms_from_ids(unique_go_b_ids),
            "shared": get_terms_from_ids(shared_go_ids)
        },
        "ppi_analysis": {
            "internal_up_regulated": internal_a[ppi_cols].to_dict('records'),
            "internal_down_regulated": internal_b[ppi_cols].to_dict('records'),
            "cross_talk": cross_talk[ppi_cols].to_dict('records')
        }
    }