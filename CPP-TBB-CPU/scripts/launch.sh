#!/bin/bash
#SBATCH --job-name=gpu_run
#SBATCH --partition=pvc-shared
#SBATCH --error=job.%J.err 
#SBATCH --output=job.%J.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=asenjo@uma.es

source /opt/intel/oneapi/setvars.sh
python ./run_minrad_levelIDC.py

