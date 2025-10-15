import pandas as pd
import os
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from contextlib import asynccontextmanager
from sqlalchemy import create_engine
from supabase import create_client, Client
import re
# **NEW:** Import the CORSMiddleware
from fastapi.middleware.cors import CORSMiddleware

# --- SUPABASE CONFIGURATION ---
# These should be set as Environment Variables in Render, not hardcoded.
SUPABASE_URL = os.environ.get("SUPABASE_URL", "https://cvtnomojgmdufcahnmkp.supabase.co")
SUPABASE_SERVICE_KEY = os.environ.get("SUPABASE_SERVICE_KEY", "YOUR_KEY_HERE")
DB_PASSWORD = os.environ.get("DB_PASSWORD", "YOUR_PASSWORD_HERE")


# --- DATA STORAGE ---
# Global variables to hold our data in memory
gene_info_df = None
go_terms_map_df = None
ensembl_to_go_df = None

# --- SETUP: LOAD DATA FROM SUPABASE ON STARTUP ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Server starting up...")
    print("Loading data from Supabase...")

    global gene_info_df, go_terms_map_df, ensembl_to_go_df
    
    try:
        # Initialize database engine
        db_url = f"postgresql://postgres:{DB_PASSWORD}@db.{SUPABASE_URL.split('//')[1]}:5432/postgres"
        engine = create_engine(db_url)
        
        # Load data from Supabase tables
        gene_info_df = pd.read_sql_table('genes', engine)
        go_terms_map_df = pd.read_sql_table('go_terms', engine)
        ensembl_to_go_df = pd.read_sql_table('ensembl_to_go', engine)
        
        # Data cleaning and uppercase conversion for consistent matching
        if 'gene_symbol' in gene_info_df.columns:
            gene_info_df['gene_symbol'] = gene_info_df['gene_symbol'].str.upper()
        
        print("✅ Data loading from Supabase complete.")
        print(f"   - Genes: {len(gene_info_df)} rows")
        print(f"   - GO Terms: {len(go_terms_map_df)} rows") 
        print(f"   - Gene-GO Mappings: {len(ensembl_to_go_df)} rows")
        
    except Exception as e:
        print(f"❌ Error loading data from Supabase: {e}")
        # In a production app, you might want to handle this more gracefully
        # For now, we let it fail loudly so we see the error in the logs.
        raise e

    print("Server is ready.")
    yield
    print("Server shutting down...")

# --- INITIALIZE FASTAPI APP ---
app = FastAPI(lifespan=lifespan)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"], # Allows all origins
    allow_credentials=True,
    allow_methods=["*"], # Allows all methods
    allow_headers=["*"], # Allows all headers
)

# --- DEFINE INPUT MODEL ---
class GeneLists(BaseModel):
    up_regulated: list[str]
    down_regulated: list[str]

# --- HEALTH CHECK ENDPOINT ---
@app.get("/")
async def root():
    return {"message": "Gene Comparison API is running", "status": "healthy"}

@app.get("/health")
async def health_check():
    data_status = {
        "genes_loaded": gene_info_df is not None and len(gene_info_df) > 0,
        "go_terms_loaded": go_terms_map_df is not None and len(go_terms_map_df) > 0,
        "mappings_loaded": ensembl_to_go_df is not None and len(ensembl_to_go_df) > 0
    }
    return {"status": "healthy", "data_loaded": data_status}

# --- API ENDPOINT ---
@app.post("/compare")
async def compare_gene_lists(lists: GeneLists):
    # Check if data is loaded
    if gene_info_df is None or go_terms_map_df is None or ensembl_to_go_df is None:
        raise HTTPException(status_code=503, detail="Service is starting up or data failed to load, please try again in a moment.")
    
    # 1. Clean and validate the input Ensembl ID lists
    ensembl_ids_a = {s.strip().upper() for s in lists.up_regulated if s.strip()}
    ensembl_ids_b = {s.strip().upper() for s in lists.down_regulated if s.strip()}

    if not ensembl_ids_a and not ensembl_ids_b:
        raise HTTPException(status_code=400, detail="Both gene lists are empty")

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

    # 4. Get gene symbols for potential future use
    symbols_a = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_a)]['gene_symbol'])
    symbols_b = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_b)]['gene_symbol'])

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
        "gene_symbols": {
            "up_regulated": list(symbols_a),
            "down_regulated": list(symbols_b)
        }
    }
