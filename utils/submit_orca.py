### IMPORTS
import sqlite3
from pathlib import Path
from extract_energies import extract_energy
import argparse
import subprocess
import threading
import time
### END

### FUNCTIONS
def is_job_running(slurm_job_id):
    result = subprocess.run(["sacct", "-j", str(slurm_job_id), "--format=State", "-X", "-n"], capture_output=True, text=True)
    if any(state in result.stdout for state in ["PREEMPTED", "DEADLINE", "FAILED", "CANCELLED", "BOOT_FAIL", "NODE_FAIL", "OUT_OF_MEMORY", "SUSPEN
DED", "TIMEOUT"]):
        return "FAIL"
    elif "PENDING" in result.stdout:
        return "PENDING"
    elif "RUNNING" in result.stdout:
        return "IN_PROCESS"
    elif "COMPLETED" in result.stdout:
        return "COMPLETED"
    else:
        return "FAIL"

def monitor_orca_jobs(conn, cur, job_dir):
    cur.execute("""
                SELECT id, inp_file, slurm_job_id
                FROM jobs
                WHERE status='orca_submitted'
                """)
    rows = cur.fetchall()
    if not rows:
        print("No submitted orca tasks yet, take a nap!!", flush=True)
        return
    for job_id, inp_file, slurm_job_id in rows:
        inp_path = Path(inp_file)
        out_file = job_dir / inp_path.name.replace(".inp",".out")
        err_file = job_dir / inp_path.name.replace(".inp",".err")

        state = is_job_running(slurm_job_id)
       
       if state == 'PENDING' or state == 'IN_PROCESS':
            continue
       
       if state == 'FAIL':
            cur.execute("""
                        UPDATE jobs
                        SET status='error', error='Slurm Error'
                        where id = ?
                        """, (job_id,))
            conn.commit()
            continue
        
        if state == "COMPLETED":
            if not out_file.exists():
                print(f"waiting for {out_file}", flush=True)
                continue
            out_text = out_file.read_text(error="ignore")

            if "ORCA TERMINATED NORMALLY" not in out_text:
                if "Error TMatrixContainers::AddMatrix" in out_text:
                    cur.execute("""
                                UPDATE jobs SET status='error', error='Error TMatrixContainers (Memory)'
                                WHERE id = ?
                                """, (job_id,))
                    conn.commit()
                    continue
                elif "multiplicity (1) is odd" in out_text:
                    cur.execute("""
                                UPDATE jobs SET status='error', error='Multiplicity error'
                                WHERE id= ?
                                """, (job_id,))
                    conn.commit()
                    continue
                else:
                    print(f"Orca failed :-( {inp_file}")
                    cur.execute("""
                            UPDATE jobs
                            SET status='error', error='orca_error'
                            WHERE id = ?
                            """, (job_id,))
                    conn.commit()
                    continue
            
            try: #If no error and job complete
                print("Trying to extract energy", flush=True)
                df = extract_energy(out_file)
                print(df, flush=True)
                
                try:
                    gibbs_free = df['Final Gibbs free energy'].iloc[0]
                except:
                    gibbs_free = None
                try:
                    sp_energy = df['FINAL SINGLE POINT ENERGY'].iloc[0] 
                except:
                    sp_energy = None
                print("gibbs", gibbs_free, "sp", sp_energy, flush=True)
                cur.execute("PRAGMA table_info(jobs)")
                columns = [info[1] for info in cur.fetchall()]
                if 'FreeEnergy' not in columns: #TODO: REMOVE THIS CHECK CAN CREATE RACE CONDITION
                    cur.execute("ALTER TABLE jobs ADD COLUMN FreeEnergy REAL")
                if 'SINGLEPOINT' not in columns:
                    cur.execute("ALTER TABLE jobs ADD COLUMN SINGLEPOINT REAL")
                conn.commit()
                cur.execute("""
                            UPDATE jobs
                            SET FreeEnergy = ?, SINGLEPOINT = ?,
                            status="orca_completed"
                            WHERE inp_file = ?
                            """, (gibbs_free, sp_energy, inp_file))
                conn.commit()
            
            except Exception as e:
                cur.execute("""
                            UPDATE jobs
                            SET status="orca_read_failed"
                            WHERE inp_file = ?
                            """, (inp_file,))
                conn.commit()
                print("Error in parsing", e)


def submit_orca(inp_file, job_dir, cpus = 4, mem = 8000):
    inp_stem = Path(inp_file).stem
    inp_file = Path(inp_file)
    slurm_filename = job_dir / f"submit_{inp_stem}.slurm"
    sbatch_script = f"""#!/bin/bash
#SBATCH --job-name={Path(inp_file).stem}
#SBATCH --nodes=1
#SBATCH --ntasks={cpus}
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time=24:00:00
#SBATCH --partition=compchem
##SBATCH --nodelist=node236,node237,node238,node239
#SBATCH --no-requeue

#Move to scratch dir for calculations

exec > >(tee -a $SCRATCH/{inp_file.name.replace(".inp", ".out")})
exec 2> >(tee -a $SCRATCH/{inp_file.name.replace(".inp", ".err")})

export PATH=/groups/kemi/julius/orca_6_1_0:$PATH
export LD_LIBRARY_PATH=/groups/julius/orca_6_1_0:$LD_LIBRARY_PATH

export PATH=/software/kemi/openmpi/openmpi-4.1.1/bin:$PATH
export LD_LIBRARY_PATH=/software/kemi/openmpi/openmpi-4.1.1/lib:$LD_LIBRARY_PATH
#module load mpi/openmpi-x86_64
cp $SLURM_SUBMIT_DIR/{inp_file} $SCRATCH
cd $SCRATCH

#/groups/kemi/hteahan/opt/orca_6_1_0_linux_x86-64_shared_openmpi418/orca {inp_file.name} --use-hwthread-cpus
/groups/kemi/julius/orca_6_1_0/orca {inp_file.name} --use-hwthread-cpus

cp *.out *.err *.xyz *trj* {job_dir}
""".format(mem=mem, cpus=cpus, job_dir=job_dir, inp_file=inp_file)
    slurm_filename.write_text(sbatch_script)
    result = subprocess.run(["sbatch", str(slurm_filename)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to submit orca job: {result.stderr}")
    job_id = result.stdout.strip().split()[-1]
    return job_id # Track slurm task by id


def monitoring_thread(DB_PATH, job_dir, poll_interval=30):
    thread_conn = sqlite3.connect(DB_PATH, timeout=60)
    thread_conn.execute("PRAGMA journal_mode=WAL;")
    thread_conn.execute("PRAGMA synchronous=NORMAL;")
    thread_cur = thread_conn.cursor()
    while True: # Initially check whether anything relevant exists - either not running yet, or missing completion.
        thread_cur.execute("""
        SELECT COUNT(*)
        FROM jobs
        WHERE (inp_file IS NOT NULL AND job_type='OPT' AND (status='inputs_created' OR status='orca_submitted')) 
        """)
        pending = thread_cur.fetchone()[0]

        if pending == 0:
            print("No more jobs pending, closing ...", flush=True)
            break
        
        monitor_orca_jobs(thread_conn, thread_cur, job_dir)  # the function from before
        time.sleep(poll_interval)


def main(args):
    
    DB_PATH = Path(args.DB_PATH)
    INP_DIR = Path(args.INP_DIR)
    OUT_DIR = Path(args.OUT_DIR)
    INP_DIR.mkdir(parents=True, exist_ok=True)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH, timeout=30)
    cur = conn.cursor()
    conn.execute("PRAGMA journal_mode=WAL;")
    conn.execute("PRAGMA synchronous=NORMAL;")

    print(f"Connection to DB made {DB_PATH}, {INP_DIR}", flush=True)
    
    t = threading.Thread(target=monitoring_thread, args=(DB_PATH,OUT_DIR))
    t.daemon = False
    t.start()
    
    while True:
    
        cur.execute("""
        SELECT COUNT(*)
        FROM jobs
        WHERE status = 'orca_submitted'
        """)
        
        running_tasks = cur.fetchone()[0]
        
        if running_tasks > args.num_parallel_tasks: # Limited to using 30 x 8 cpus so I dont take up too much.
            print("Too many tasks running, sleeping...", running_tasks, flush=True)
            time.sleep(30)
            continue
        
        # Fetch one task
        cur.execute("""
        UPDATE jobs
        SET status='orca_claimed'
        WHERE id = (
            SELECT id
            FROM jobs
            WHERE job_type='OPT'
            AND status='inputs_created'
            LIMIT 1)
        RETURNING id, inp_file, ligand_id;
        """)

        row = cur.fetchone()
        
        if row is None:
            time.sleep(5)
            continue

        job_id, inp_file, ligand_id = row
        conn.commit()
               
        try:
            OPT_inp_file = str(inp_file).replace("INPUT_SP", "INPUT")
            job_id = submit_orca(OPT_inp_file, job_dir=OUT_DIR, cpus=args.cpus, mem=args.mem)
        
            print(f"Submitted ORCA job {OPT_inp_file} as Slurm ID {job_id}", flush=True)
            
            cur.execute("""
            UPDATE jobs
            SET status='orca_submitted', slurm_job_id=?, inp_file=?
            WHERE inp_file=?
            """, (job_id, str(OPT_inp_file), str(inp_file)))
            
            conn.commit()
        
        except Exception as e:
            
            print(f"Error submitting {inp_file}: {e}")
            
            cur.execute("""
                UPDATE jobs
                SET status='error_in_orca_submit', error=?
                WHERE inp_file=?
            """, (str(e), str(inp_file)))
        
            conn.commit()
        
            continue
    t.join()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--DB_PATH", type=str, help="Path to DB file name. Input example: .../Projects/db/jobs.db. PATH MUST ALREADY EXIST (created by seed_db.py)", required=True)
    parser.add_argument("--INP_DIR", type=str, help="Path to Input file directory", required=True)
    parser.add_argument("--OUT_DIR", type=str, help="Path to Output file directory", required=True)
    parser.add_argument("--cpus", type=int, help="CPUs available for ORCA calculation (Make sure it is consistent with input file)", default=8)
    parser.add_argument("--mem", type=int, help="Memory available for ORCA calculation", default=20000)
    parser.add_argument("--num_parallel_tasks", type=int, help="Number of parallel calculations to run", default=25)
    args = parser.parse_args()
    print("Beginning", flush=True)
    
    main(args)
        





























