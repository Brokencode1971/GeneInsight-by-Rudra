import pandas as pd
import os
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from contextlib import asynccontextmanager
from sqlalchemy import create_engine
from supabase import create_client, Client
import re
# **THIS IS THE FIX:** Import the CORSMiddleware
from fastapi.middleware.cors import CORSMiddleware

# --- SUPABASE CONFIGURATION ---
SUPABASE_URL = os.environ.get("SUPABASE_URL")
SUPABASE_SERVICE_KEY = os.environ.get("SUPABASE_SERVICE_KEY")
DB_PASSWORD = os.environ.get("DB_PASSWORD")

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
        
    except Exception as e:
        print(f"❌ Error loading data from Supabase: {e}")
        raise e

    print("Server is ready.")
    yield
    print("Server shutting down...")

# --- INITIALIZE FASTAPI APP ---
app = FastAPI(lifespan=lifespan)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# --- DEFINE INPUT MODEL ---
class GeneLists(BaseModel):
    up_regulated: list[str]
    down_regulated: list[str]

# --- HEALTH CHECK ENDPOINT ---
@app.get("/")
async def root():
    return {"message": "Gene Comparison API is running", "status": "healthy"}

# --- API ENDPOINT ---
@app.post("/compare")
async def compare_gene_lists(lists: GeneLists):
    if gene_info_df is None or go_terms_map_df is None or ensembl_to_go_df is None:
        raise HTTPException(status_code=503, detail="Service starting up, please try again in a moment")
    
    ensembl_ids_a = {s.strip().upper() for s in lists.up_regulated if s.strip()}
    ensembl_ids_b = {s.strip().upper() for s in lists.down_regulated if s.strip()}

    if not ensembl_ids_a and not ensembl_ids_b:
        raise HTTPException(status_code=400, detail="Both gene lists are empty")

    db_ensembl_ids = set(gene_info_df['ensembl_gene_id'].dropna())
    mapped_a_count = len(ensembl_ids_a.intersection(db_ensembl_ids))
    mapped_b_count = len(ensembl_ids_b.intersection(db_ensembl_ids))

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

    symbols_a = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_a)]['gene_symbol'])
    symbols_b = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_b)]['gene_symbol'])

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
