import pandas as pd
import os
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from contextlib import asynccontextmanager
from sqlalchemy import create_engine
from fastapi.middleware.cors import CORSMiddleware
import re

# --- SUPABASE CONFIGURATION (from Environment Variables) ---
SUPABASE_URL = os.environ.get("SUPABASE_URL")
DB_PASSWORD = os.environ.get("DB_PASSWORD")

# --- DATA STORAGE & DB CONNECTION ---
# Global variables to hold the database engine
engine = None

# --- SETUP: CREATE DATABASE CONNECTION ON STARTUP ---
@asynccontextmanager
async def lifespan(app: FastAPI):
    print("Server starting up...")
    print("Establishing database connection...")
    global engine
    try:
        if not SUPABASE_URL or not DB_PASSWORD:
            raise ValueError("SUPABASE_URL and DB_PASSWORD environment variables must be set.")
            
        db_url = f"postgresql://postgres:{DB_PASSWORD}@db.{SUPABASE_URL.split('//')[1]}:5432/postgres"
        engine = create_engine(db_url)
        # Test the connection
        with engine.connect() as connection:
            print("✅ Database connection successful.")
    except Exception as e:
        print(f"❌ FATAL: Could not connect to the database: {e}")
        engine = None # Ensure engine is None if connection fails
    
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
    if engine is None:
        raise HTTPException(status_code=503, detail="Database connection is not available.")

    with engine.connect() as connection:
        # 1. Clean and validate the input Ensembl ID lists
        ensembl_ids_a = {s.strip().upper() for s in lists.up_regulated if s.strip()}
        ensembl_ids_b = {s.strip().upper() for s in lists.down_regulated if s.strip()}

        # 2. Get the count of mapped genes
        gene_info_df = pd.read_sql("SELECT ensembl_gene_id, gene_symbol FROM genes", connection)
        db_ensembl_ids = set(gene_info_df['ensembl_gene_id'].str.upper().dropna())
        mapped_a_count = len(ensembl_ids_a.intersection(db_ensembl_ids))
        mapped_b_count = len(ensembl_ids_b.intersection(db_ensembl_ids))

        # 3. Perform GO Term Comparison
        ensembl_to_go_df = pd.read_sql("SELECT ensembl_gene_id, go_id FROM ensembl_to_go", connection)
        go_a = ensembl_to_go_df[ensembl_to_go_df['ensembl_gene_id'].isin(ensembl_ids_a)]
        go_b = ensembl_to_go_df[ensembl_to_go_df['ensembl_gene_id'].isin(ensembl_ids_b)]
        
        go_set_a = set(go_a['go_id'])
        go_set_b = set(go_b['go_id'])
        
        unique_go_a_ids = list(go_set_a - go_set_b)
        unique_go_b_ids = list(go_set_b - go_set_a)
        shared_go_ids = list(go_set_a & go_set_b)

        go_terms_map_df = pd.read_sql("SELECT \"GO_ID\", \"GO_Term\" FROM go_terms", connection)
        def get_terms_from_ids(id_list):
            if not id_list: return []
            terms = go_terms_map_df[go_terms_map_df['GO_ID'].isin(id_list)]
            return terms.rename(columns={'GO_ID': 'id', 'GO_Term': 'term'}).to_dict('records')

        # 4. Perform PPI Network Analysis
        symbols_a = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_a)]['gene_symbol'].str.upper())
        symbols_b = set(gene_info_df[gene_info_df['ensembl_gene_id'].isin(ensembl_ids_b)]['gene_symbol'].str.upper())
        
        biogrid_ppi_df = pd.read_sql("SELECT \"Official Symbol Interactor A\", \"Official Symbol Interactor B\" FROM biogrid_ppi", connection)
        biogrid_ppi_df['Official Symbol Interactor A'] = biogrid_ppi_df['Official Symbol Interactor A'].str.upper()
        biogrid_ppi_df['Official Symbol Interactor B'] = biogrid_ppi_df['Official Symbol Interactor B'].str.upper()

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
