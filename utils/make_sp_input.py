### IMPORTS
from xyz2inp import xyz_to_inp
import argparse
import sqlite3
import json
from pathlib import Path
### END


### FUNCTIONS

def main(args):
    DB_PATH = Path(args.DB_PATH)
    INP_DIR = Path(args.INP_DIR)
    INP_DIR.mkdir(parents=True, exist_ok=True)
    EMBED_DIR = Path(args.EMBED_DIR)
    conn = sqlite3.connect(DB_PATH, timeout=30)
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA busy_timeout = 30000;")
    cur = conn.cursor()
    print(f"Connection to DB made {DB_PATH}, {INP_DIR}", flush=True)
     
    while True:
        xyz_file, ligand_id = None, None
        
        #Atomic selection of embedded molecules
        try:
            cur.execute("""
            UPDATE jobs
            SET status='inputting'
            WHERE id = (
                SELECT id
                FROM jobs
                WHERE job_type='embedding'
                AND status='complete'
                LIMIT 1)
            RETURNING id, xyz_file, ligand_id, xTBEnergy;
            """)

            row = cur.fetchone()
            
            if row is None:
                # Free the reserved input
                
                time.sleep(5)
                continue
            
            conn.commit()

            job_id, xyz_file, ligand_id, xTB_E= row
            
            # xyz files and energies stored as json
            xyz_file = json.loads(xyz_file)
            xTB_E = json.loads(xTB_E)
            
            charge, multiplicity = args.charge, args.multiplicity
        
            for i, xyz in enumerate(xyz_file):
                xyzsfilepath = EMBED_DIR / xyz
                inp_file = xyz_to_inp(xyzsfilepath, charge, multiplicity, input_name=INP_DIR, keep_lines=2, query = f"{args.query} \n%pal \nnprocs {args.cpus} \nend \n%maxcore {int(args.mem/args.cpus*0.75)}", please_return=True)
            
                if not inp_file:
                    raise RuntimeError("No ORCA input files generated")
                
                # Create new row
                cur.execute("""
                            INSERT INTO jobs (job_type, xyz_file, xTBEnergy, ligand_id, inp_file, status)
                            VALUES (?,?,?,?,?,?)
                            """, ('SP', str(xyz), xTB_E[i], ligand_id, str(inp_file), 'pending'))
            
            # Update the old selection
            cur.execute("""
                UPDATE jobs
                SET status='inputs_created'
                WHERE id=?
                """, (job_id,))
            conn.commit()
            continue

        except Exception as e:

            conn.rollback()
            
            print("Error in INP FILE generation:", e, flush=True)
            
            if xyz_file is not None and ligand_id is not None:
                try:
                    cur.execute("""
                    UPDATE jobs
                    SET status='input_gen_failed'
                    WHERE id=?
                    """, (job_id,))
                    conn.commit()
                except:
                    conn.rollback()
            continue   
        
        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--DB_PATH", type=str, help="Path to DB file name. Input example: .../Projects/db/jobs.db. PATH MUST ALREADY EXIST (created by seed_db.py)", required=True)
    parser.add_argument("--EMBED_DIR", type=str, help="Path to .xyzs file directory from embedding", required=True)
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
