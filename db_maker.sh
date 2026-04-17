#!/bin/bash
micromamba run -n Complex python ../../scripts/db_seed.py --SMILES Missing_CSD_N.csv --DB_PATH db/jobs.db
