### IMPORTS
from xyz2inp import xyz_to_inp
import argparse
import sqlite3
import json
from pathlib import Path
import multiprocessing
from multiprocessing import Pool
import os
import time
### END


### FUNCTIONS
def select_create(conn, args):
    cur = conn.cursor()
    cur.execute("""BEGIN IMMEDIATE""")
    
    cur.execute("""
                UPDATE jobs
                SET status='inputting'

                WHERE id = (
                    SELECT id
                    FROM jobs
                    WHERE job_type='SP'
                    AND status='orca_completed'
                    AND SINGLEPOINT is not null
                    AND ligand_id NOT in (
                        SELECT ligand_id FROM jobs WHERE status='inputting'
                )
                    order by ligand_id ASC, SINGLEPOINT ASC, id ASC
                    LIMIT 1)
                RETURNING id, xyz_file, ligand_id, SINGLEPOINT
                """)
    
    row = cur.fetchone()
    
    if row is None:
        time.sleep(5)
        conn.rollback()
        return

    conn.commit()
    
    try:
    
        job_id, xyz_file, ligand_id, SP_E = row
        xyzsfilepath = Path(args.EMBED_DIR) / xyz_file

        charge, multiplicity = args.charge, args.multiplicity
    
    
        inp_file = xyz_to_inp(str(xyzsfilepath),
                                    charge,
                                    multiplicity,
                                    input_name=str(args.INP_DIR),
                                    keep_lines=2,
                                    query=f"{args.query} \n %pal \nprocs {args.cpus} \nend \n%maxcore {int(args.mem/args.cpus*0.75)}",
                                    please_return=True
                            )

        if not inp_file:
            raise RuntimeError("No ORCA input files generated")
        
        cur.execute("""
                    UPDATE jobs
                    SET status='inputs_created', job_type='OPT'
                    WHERE id=?
                    """,
                    (job_id,))
        conn.commit()
    
    except Exception as e:
    
        conn.rollback()

        print("Error in INP FILE generation", e, flush=True)
        if job_id is not None:
            try:
                cur.execute("""
                            UPDATE jobs
                            SET status='input_gen_failed', error=?
                            WHERE id=?
                            """, 
                            (str(e), job_id,))
                conn.commit()
            except Exception as e:
                print("FAILED TO MARK JOB AS FAILED", e, flush=True)
                return
     

def worker(args):
    conn = sqlite3.connect(args.DB_PATH, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA busy_timeout = 30000;")
    
    while True:
        select_create(conn, args)
    


def main(args):
    DB_PATH = Path(args.DB_PATH)
    INP_DIR = Path(args.INP_DIR)
    INP_DIR.mkdir(parents=True, exist_ok=True)
    
    n_cpus = int(os.environ.get("SLURM_CPUS_PER_TASK", 1))
    print("Running using", n_cpus, " CPUS")
    pool = Pool(processes=n_cpus)
    
    with pool as p:
        _ = p.map(worker, [args])
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--DB_PATH", type=str, help="Path to DB file name. Input example: .../Projects/db/jobs.db. PATH MUST ALREADY EXIST (created by seed_db.py)", required=True)
    parser.add_argument("--EMBED_DIR", type=str, help="Path to .xyzs file directory from embedding", default = "/groups/kemi/hteahan/Orca_workflow/Automated/Pipeline/Embed")
    parser.add_argument("--INP_DIR", type=str, help="Path to Input file directory", required=True)
    parser.add_argument("--query", type=str, help="First line query in orca inp file", default = "!r2scan-3c freq opt")
    parser.add_argument("--charge", type=int, help="Charge on complex", default=0)
    parser.add_argument("--multiplicity", type=int, help="Multiplicity of complex", default=1)
    parser.add_argument("--cpus", type=int, help="CPUs available for ORCA calculation (Written in input file)", default=8)
    parser.add_argument("--mem", type=int, help="Memory avaible for ORCA calculation (Written in input file)", default=20000)
    parser.add_argument("--override", type=str, help="Path to text file with ligand_id's to override", default=None)
    args = parser.parse_args()
    print("Beginning", flush=True)
    main(args)
