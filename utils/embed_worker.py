### IMPORTS
import sqlite3
from pathlib import Path
import json
from embed import embed
import time
import argparse
### END

### FUNCTIONS
def main(args):
    DB_PATH = Path(args.DB_PATH)
    EMBED_DIR = Path(args.EMBED_DIR)
    EMBED_DIR.mkdir(parents=True, exist_ok=True) # Directory where embedding files go?
    
    conn = sqlite3.connect(DB_PATH, timeout=30)
    #conn.execute("PRAGMA journal_mode=WAL;")
    #conn.execute("PRAGMA synchronous=NORMAL;")

    cur = conn.cursor()
    print(f"Connection to DB made {DB_PATH}, {EMBED_DIR}", flush=True)
    #cur.execute("""
    #BEGIN IMMEDIATE;
    #""")
    
    while True:

        cur.execute("""
        UPDATE jobs
        SET status='running'
        WHERE id = (
            SELECT id
            FROM jobs
            WHERE job_type='embedding'
            AND status='pending'
            LIMIT 1)
        RETURNING id, smiles, ligand_id;
        """)

        row = cur.fetchone()
        print(row)
        if row is None:
            conn.rollback() # No more pending tasks
            print("No pending embedding jobs.", flush=True)
            return

        job_id, smile, ligand_id = row
        conn.commit()


        cur.execute("""
        UPDATE jobs
        SET status='running'
        WHERE id=?
        """, (job_id,))

        conn.commit()
        
        try:
            print(smile, job_id, ligand_id, args.xTB_path, EMBED_DIR, args.N_tries, args.cpus)
            conformers = embed(smile, job_id, ligand_id, args.xTB_path, EMBED_DIR, args.N_tries, cpus = args.cpus)
            print(f"embedded {job_id} : {ligand_id}", conformers, flush=True)
            cur.execute("""
            UPDATE jobs
            SET status='complete',
            xTBEnergy=?,
            xyz_file=?
            WHERE id=?
            """, (json.dumps(conformers["xTB_energies"]),json.dumps(conformers["XYZ"]), job_id))
            conn.commit()
        except Exception as e:
            print(f"Exception in {job_id} : {ligand_id}: {e}", flush=True)
            cur.execute("""
            UPDATE jobs
            SET status='failed',
            error=?
            WHERE id=?
            """, (str(e), job_id))
            conn.commit()

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--DB_PATH", type=str, help="Path to DB file name. Input example: .../Projects/db/jobs.db. PATH MUST ALREADY EXIST (created by seed_db.py)", required=True)
    parser.add_argument("--EMBED_DIR", type=str, help="Path to Embedding directory", required=True)
    parser.add_argument("--xTB_path", type=str, help="Full path to xTB binary", default="/groups/kemi/hteahan/opt/xtb-6.7.1/xtb-dist/bin/xtb")
    parser.add_argument("--N_tries", type=int, help="Number of embedding attempts", default=1)
    parser.add_argument("--cpus", type=int, help="Number of CPUs", default=2)
    args = parser.parse_args()
    print("Beginning", flush=True)

    main(args)
