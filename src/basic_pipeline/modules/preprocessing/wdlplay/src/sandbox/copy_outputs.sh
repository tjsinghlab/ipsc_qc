#!/bin/bash
#SBATCH --job-name=copy_caper_outputs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=100:00:00
#SBATCH --output=output.txt
#SBATCH --error=error.txt

source activate wdlplay

# start=$(( ($SLURM_ARRAY_TASK_ID - 1) * 100 ))
# end=$(( $SLURM_ARRAY_TASK_ID * 100 ))

python /gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/src/sandbox/copy_outputs.py