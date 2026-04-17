### IMPORTS
import sqlite3
import pandas as pd
import argparse 
import json
import time
from pathlib import Path
### END

### FUNCTIONS:
def main(args):
    """ INITIALIZES SQLITE3 DB
        - Loads csv containing smiles & information regarding complexes
        - Creates Sqlite3 DB with pending status
    """
    DB_PATH = Path(args.DB_PATH)
    DB_PATH.parent.mkdir(exist_ok=True, parents=True)

    print(f"Making database at {DB_PATH}", flush=True)
    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;") # NOTE: Allows for writing from both monitoring thread & main concurrently using Wriet-Ahead-Logging
    cur = conn.cursor()

    # Create fields:
    cur.execute("""
    CREATE TABLE IF NOT EXISTS jobs (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_type TEXT NOT NULL,
    ligand_id INTEGER,
    mol_id INTEGER,
    smiles TEXT,
    ligand_smiles TEXT,
    xyz_file TEXT,
    inp_file TEXT,
    out_file TEXT,
    xTBEnergy REAL,
    FreeEnergy REAL,
    status TEXT NOT NULL,
    slurm_job_id TEXT,
    error TEXT,

    created_at DATETIME DEFAULT CURRENT_TIMESTAMP,
    updated_at DATETIME DEFAULT CURRENT_TIMESTAMP
    );
    """)    
    SMILES_CSV = args.SMILES
    df = pd.read_csv(SMILES_CSV)

    if 'csd_ID' in df.columns():
        for _, row in df.iterrows():
            cur.execute("""
                INSERT INTO jobs (job_type, ligand_id, csd_ID, smiles, ligand_smiles, status)
                VALUES ('embedding', ?, ?, ?, ?, 'pending')
            """, (row["ligand_ID"], row['csd_ID'], row[args.SMILES_COL], row[args.LIGAND_SMILES_COL]))
    else:
        for _, row in df.iterrows():
            cur.execute("""
                INSERT INTO jobs (job_type, ligand_id, smiles, ligand_smiles, status)
                VALUES ('embedding', ?, ?, ?, 'pending')
            """, (row["ligand_ID"], row[args.SMILES_COL], row[args.LIGAND_SMILES_COL]))
 
    conn.commit()
    conn.close()

    print("DB seeded with embedding jobs.", flush=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--SMILES", type=str, help=".csv path containing smiles. Must contain ligand_id, mol_id.", required=True)
    parser.add_argument("--SMILES_COL", type=str, help="Column name for column containing the smiles of interest.", default="smiles")
    parser.add_argument("--DB_PATH", type=str, help="Path to DB file name. Input example: .../Projects/db/jobs.db (Both db and jobs.db will be created)", default="db/jobs.db")
    parser.add_argument("--LIGAND_SMILES_COL", type=str, help="Column name for column containing the ligands corresponding to the complexes.", default="ligand_smiles")
    args = parser.parse_args()
    try:
        main(args)
    except Exception as e:
        print(e, flush=True)

