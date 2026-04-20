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
    result = subprocess.run(["sacct", "-j", str(slurm_job_id), "--format=JobID,State","--array",  "-X", "-n"], capture_output=True, text=True)
    
    lines = result.stdout.strip().splitlines()

    states = {}

    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        job_id, state = parts[0], parts[1]
        #print("job_id", job_id, "state", state, flush=True)
        if job_id.endswith((".batch", ".extern")):
            continue
        if any(s in state for s in ["PREEMPTED", "DEADLINE", "FAILED", "CANCELLED", "BOOT_FAIL", "NODE_FAIL", "OUT_OF_MEMORY", "SUSPENDED", "TIMEOUT"]):
            states[job_id] = "FAIL"
            print("FAIL", state, flush=True)
        elif state == "RUNNING":
            states[job_id] = "IN_PROCESS"
        elif state == "PENDING":
            states[job_id] = "PENDING"
        elif state == "COMPLETED":
            states[job_id] = "COMPLETED"
#            print("COMPLETED", state, flush=True)
        else:
            states[job_id] = "FAIL"
    return states
        
def monitor_orca_jobs(conn, cur, job_dir, columns, tracker):

    cur.execute("""
                SELECT id, inp_file, slurm_job_id, slurm_task_id, error
                FROM jobs
                WHERE status='orca_sp_submitted'
                """)
    rows = cur.fetchall()
    print("Number of submitted tasks:", len(rows), flush=True) 
    if not rows:
        print("No submitted sp orca tasks yet, take a nap!!", flush=True)
        return

    grouped = {}

    for row in rows:
        db_id, inp_file, slurm_job_id, task_id, error = row

        if slurm_job_id not in grouped:
            grouped[slurm_job_id] = []

        grouped[slurm_job_id].append(row)
    try:        
        for slurm_job_id, group_rows in grouped.items():
            
            states = is_job_running(slurm_job_id)
            for db_id, inp_file, slurm_job_id, task_id, error in group_rows:
                inp_path = Path(inp_file)
                out_file = job_dir / inp_path.name.replace(".inp",".out")
                err_file = job_dir / inp_path.name.replace(".inp",".err")
                full_task_id = f"{slurm_job_id}_{task_id}"
    #            if task_id=="":
     #               state = states.get(slurm_job_id, "UNKNOWN")
     #           else:
                
                state = states.get(full_task_id, "UNKNOWN")
                print(full_task_id, state, flush=True)

                if state == "UNKNOWN":
                    new_state = states.get(slurm_job_id, "UNKNOWN")
                    if new_state != "UNKNOWN" and new_state != None:
                        state = new_state
                        print("Used new state",db_id,slurm_job_id,task_id, flush=True)

                #print("job_id", db_id, "state", state, flush=True)
                if state == 'PENDING' or state == 'IN_PROCESS':
                    print("state pending", db_id, state, inp_file, flush=True)
                    continue 
                if state == 'FAIL':
                    print("Failed", db_id, flush=True)
                    cur.execute("""
                                UPDATE jobs
                                SET status='error', error='Slurm Error'
                                where id = ?
                                """, (db_id,))
                    conn.commit()
                    
                    
                if state == "COMPLETED":
                    print("Completed - ", db_id,slurm_job_id, task_id, flush=True)
                    if not out_file.exists():
                        tracker[db_id] = tracker.get(db_id, 0) + 1
                        print("n_attempts", tracker.get(db_id, 0)+1)
                        if tracker[db_id] > 10:
                            print("Couldn't find output file", flush=True)
                            cur.execute("""
                                        UPDATE jobs
                                        SET status='error', error='Output_not_found'
                                        WHERE id=?
                                        """, (db_id,))
                            conn.commit()
                            continue
                        else:
                            print("Retrying", flush=True)
                        continue
                        
                    out_text = out_file.read_text(errors="ignore")
                    if "ORCA TERMINATED NORMALLY" not in out_text:
                        if "Error TMatrixContainers::AddMatrix" in out_text:
                            cur.execute("""
                                        UPDATE jobs SET status='error', error='Error TMatrixContainers (Memory)'
                                        WHERE id = ?
                                        """, (db_id,))
                            conn.commit()
                            continue
                        elif "multiplicity (1) is odd" in out_text:
                            cur.execute("""
                                        UPDATE jobs SET status='error', error='Multiplicity error'
                                        WHERE id= ?
                                        """, (db_id,))
                            conn.commit()
                            continue
                        else:
                            print(f"Orca failed :-( {inp_file}")
                            cur.execute("""
                                    UPDATE jobs
                                    SET status='error', error='orca_sp_error'
                                    WHERE id = ?
                                    """, (db_id,))
                            conn.commit()
                    try: #If no error and job complete
                        print("Trying to extract SP energy", flush=True)
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
                       

                      # inp_file_new = Path(str(inp_file).replace("INPUT_SP", "INPUT"))
                        cur.execute("""
                                    UPDATE jobs
                                    SET FreeEnergy = ?, 
                                    SINGLEPOINT = ?,
                                    status='orca_completed'
                                    WHERE id = ?
                                    """, (gibbs_free, sp_energy, db_id)) 
                        conn.commit()
                        continue        
                    except Exception as e:
                        cur.execute("""
                                    UPDATE jobs
                                    SET status="orca_sp_read_failed"
                                    WHERE id = ?
                                    """, (db_id,))
                        conn.commit()
                        print("Error in parsing", e, flush=True)
    except Exception as e:
        print("ERROR", e, flush=True)
        
def submit_orca_array(inp_files, job_dir, BATCH_SIZE, partition, cpus = 4, mem = 8000):
    
    #inp_stem = [Path(inp_file).stem for inp_file in inp_files]
    inp_files = [Path(inp_file) for inp_file in inp_files]
    timestamp = int(time.time())
    
    slurm_filename = job_dir / f"submit_{timestamp}.slurm"
    batch_file = job_dir / f"batch_array_{timestamp}.txt"

    with open(batch_file, "w") as f:
        for inp in inp_files:
            f.write(str(inp) + "\n")
    
    array_len = len(inp_files)
    max_parallel = min(BATCH_SIZE, array_len)


    sbatch_script = """#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --nodes=1
#SBATCH --ntasks={cpus}
#SBATCH --cpus-per-task=1
#SBATCH --mem={mem}
#SBATCH --time=10:00:00
#SBATCH --partition={partition}
#SBATCH --array=1-{array_len}%{max_parallel}
##SBATCH --nodelist=node236,node238,node237,node239
#SBATCH --no-requeue

# Create each array entry:
INPUT=$(sed -n "$SLURM_ARRAY_TASK_ID"p {batch_file}) 
BASENAME=$(basename $INPUT .inp)

exec > >(tee -a $SCRATCH/${{BASENAME}}.out)
exec 2> >(tee -a $SCRATCH/${{BASENAME}}.err)

#export PATH=/groups/kemi/julius/orca_6_1_0:$PATH
#export LD_LIBRARY_PATH=/groups/julius/orca_6_1_0:$LD_LIBRARY_PATH

export PATH=/software/kemi/openmpi/openmpi-4.1.1/bin:$PATH
export LD_LIBRARY_PATH=/software/kemi/openmpi/openmpi-4.1.1/lib:$LD_LIBRARY_PATH

cp $INPUT $SCRATCH
cd $SCRATCH

/groups/kemi/julius/orca_6_1_0/orca $(basename $INPUT) --use-hwthread-cpus
mkdir -p {job_dir}
cp -v *.out *.err *.xyz *trj* {job_dir}/ 2>/dev/null || true
""".format(job_name=batch_file.stem, mem=mem, cpus=cpus, job_dir=job_dir, batch_file=batch_file, array_len=array_len, max_parallel=max_parallel, partition=partition)
    slurm_filename.write_text(sbatch_script)
    result = subprocess.run(["sbatch", str(slurm_filename)],
                            capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Failed to submit orca job: {result.stderr}")
    job_id = result.stdout.strip().split()[-1]

    return job_id # Track slurm task by id

def monitoring_thread(DB_PATH, job_dir, poll_interval=10):
    thread_conn = sqlite3.connect(DB_PATH, timeout=poll_interval, check_same_thread=False)
    thread_conn.execute("PRAGMA journal_mode=WAL;")
    thread_conn.execute("PRAGMA synchronous=NORMAL;")
    thread_cur = thread_conn.cursor()
    thread_cur.execute("PRAGMA table_info(jobs)")
    columns = [info[1] for info in thread_cur.fetchall()]
    tracker = {}
  
    if 'FreeEnergy' not in columns:
        thread_cur.execute("ALTER TABLE jobs ADD COLUMN FreeEnergy REAL")
    if 'SINGLEPOINT' not in columns:
        thread_cur.execute("ALTER TABLE jobs ADD COLUMN SINGLEPOINT REAL")
    thread_conn.commit()
    
    
    while True: # Initially check whether anything relevant exists - either not running yet, or missing completion.
        print("monitoring thread running")
        print(tracker, flush=True)
        

        try:
            thread_cur.execute("""
            SELECT COUNT(*)
            FROM jobs
            WHERE (inp_file IS NOT NULL AND job_type='SP' AND (status='pending' OR status='orca_sp_submitted')) 
            """)
        
            pending = thread_cur.fetchone()[0]
         
            if pending == 0:
                print("No more jobs pending, sleeping ...", flush=True)
                time.sleep(60)
            
            monitor_orca_jobs(thread_conn, thread_cur, job_dir, columns, tracker)  # monitoring completed and failed jobs
            time.sleep(poll_interval)

        except sqlite3.OperationalError as e:
            if "database is locked" in str(e):
                print("DB locked", flush=True)    
                
                #thread_conn = sqlite3.connecti(DB_PATH, timeout=60, check_same_thread=False)
                thread_conn.rollback()
                time.sleep(5)
                continue
            else:
                print(e, flush=True)
                break
                
            
    thread_conn.close()


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
   
    cur.execute("PRAGMA table_info(jobs)")
    columns = [c[1] for c in cur.fetchall()]

    if "slurm_task_id" not in columns:
        cur.execute("ALTER TABLE jobs ADD COLUMN slurm_task_id INTEGER")
        conn.commit()

    while True:
        # Throttler - allow only certain number of parallel tasks.
        #cur.execute("""BEGIN IMMEDIATE""")
        cur.execute("""
        SELECT COUNT(*)
        FROM jobs
        WHERE status = 'orca_sp_submitted'
        """)
    
        running_tasks = cur.fetchone()[0]
        
        print("Number of running tasks", running_tasks, flush=True)
        if running_tasks > args.batch_size:
            print("Too many tasks running, sleeping...", running_tasks, flush=True)
            time.sleep(20)
            continue
        
        cur.execute("""
        SELECT COUNT(*)
        FROM jobs
        WHERE job_type='SP' and status='pending'
        """)
        
        pending = cur.fetchone()[0]

        print("Pending tasks", pending, flush=True)
        if pending == 0:
            if running_tasks == 0:
                print("No pending or running tasks, quitting...", flush=True)
                break
            print("No pending tasks, sleeping...", flush=True)
            time.sleep(20)
            continue



        # Reworked to work as array process
        BATCH_SIZE = args.batch_size # args.batch_size
        # BEGIN IMMEDIATE - TODO
        cur.execute(f"""
        UPDATE jobs
        SET status='orca_claimed'
        WHERE id IN (
            SELECT id
            FROM jobs
            WHERE job_type='SP'
            AND status='pending'
            LIMIT ?)
        RETURNING id, inp_file, ligand_id;
        """, (BATCH_SIZE,))

        rows = cur.fetchall()
        
        if not rows:
            conn.rollback()
            time.sleep(10) #Wait
            continue
        
        #for row in rows:
        #    job_id, inp_file, ligand_id = row

        conn.commit()
        
        try:
            inp_files = [row[1] for row in rows]
    
            job_id = submit_orca_array(inp_files, job_dir=OUT_DIR, BATCH_SIZE=BATCH_SIZE, partition=args.partition, cpus=args.cpus, mem=args.mem)
            
            print(f"Submitted ORCA job {inp_files[0]}-{inp_files[-1]} as Slurm ID {job_id}", flush=True)
            
            for i, (db_id, inp_file, ligand_id) in enumerate(rows, start=1): 
                cur.execute("""
                UPDATE jobs
                SET status='orca_sp_submitted', slurm_job_id=?, slurm_task_id=?
                WHERE id=?
                """, (job_id, i, int(db_id)))
            conn.commit()
        
        except Exception as e:
            
            print(f"Error submitting {inp_files[0]}-{inp_files[-1]}: {e}")
            
            for db_id, inp_file, ligand_id in rows:
                    
                cur.execute("""
                    UPDATE jobs
                    SET status='error_in_orca_sp_submit', error=?
                    WHERE id=?
                """, (str(e), int(db_id)))

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
    parser.add_argument("--batch_size", type=int, help="Number of parallel array calculations to run", default=100)
    parser.add_argument("--partition", type=str, help="Partition", default="kemi1")
    args = parser.parse_args()
    print("Beginning", flush=True)
    
    main(args)
        





























