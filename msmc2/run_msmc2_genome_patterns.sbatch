#!/bin/bash
#SBATCH -p general
#SBATCH -J lutzo_msmc2
#SBATCH --cpus-per-task=1
#SBATCH -t 24:00:00
#SBATCH --mem=8g
#SBATCH --array=1-132
#SBATCH -o /work/users/d/t/dturissi/lutzomyia/msmc2/logs/run_msmc2_genome_patterns_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/lutzomyia/msmc2/logs/run_msmc2_genome_patterns_%A_%a.err

ml msmc2/2.1.3


i=$((${SLURM_ARRAY_TASK_ID} - 1))


sh_file="/work/users/d/t/dturissi/lutzomyia/msmc2/run_msmc2_genome.sh"
msmc2_cmd=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${sh_file})

${msmc2_cmd} "1*2+15*1+1*2"
${msmc2_cmd} "1*4+12*1+1*4"
${msmc2_cmd} "1*4+15*1+1*4"
${msmc2_cmd} "1*4+20*1+1*4"
${msmc2_cmd} "1*4+8*1+1*4"
${msmc2_cmd} "1*6+8*1+1*4"
${msmc2_cmd} "1*6+8*1+1*6"
${msmc2_cmd} "1*6+8*2+1*6"
${msmc2_cmd} "1*8+12*1+1*8"
${msmc2_cmd} "1*8+8*1+1*8"