#!/bin/sh
#SBATCH -c 24
#SBATCH -t 01:40:00
#SBATCH -p thinnodes,cola-corta
#SBATCH --mail-type=end
#SBATCH --mail-user=felipepower7@gmail.com
#SBATCH --output slurm_V2124py24.out
srun python ex_python3.py
