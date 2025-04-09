#!/bin/bash
#SBATCH -p general
#SBATCH -J lutzo_msmc2
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=8g
#SBATCH --array=1-132
#SBATCH -o /work/users/d/t/dturissi/lutzomyia/msmc2/logs/prep_msmc2_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/lutzomyia/msmc2/logs/prep_msmc2_%A_%a.err


i=$((${SLURM_ARRAY_TASK_ID} - 1))


sh_file="/work/users/d/t/dturissi/lutzomyia/msmc2/prep_msmc2.sh"
msmc2_cmd=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${sh_file})

${msmc2_cmd}