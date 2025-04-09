#!/bin/bash
#SBATCH -p general
#sbatch -J lutzo_pi_vcftools
#SBATCH -t 24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --array=1-8
#SBATCH -o /work/users/d/t/dturissi/lutzomyia/poly_windows/logs/get_lutzo_pi_vcftools_%A_%a.out
#SBATCH -e /work/users/d/t/dturissi/lutzomyia/poly_windows/logs/get_lutzo_pi_vcftools_%A_%a.err



ml vcftools/0.1.15


#pops=('All', 'Colombia' 'Jacobina_L_S1' 'Jacobina_M_S2' 'Laphina' 'Marajo' 'Sobral_1S_A' 'Sobral_1S_B' 'Sobral_2S' 'Jacobina' 'Sobral_1S')
pops=('All')

i=$((${SLURM_ARRAY_TASK_ID} - 1))
pop=${pops[${i}]}

#define dirs and files
base_dir=/work/users/d/t/dturissi/lutzomyia/poly_windows
output_dir=${base_dir}/pi_vcftools_results
pop_id_file=${output_dir}/${pop}_ids.txt



#make tmp dir and change to it
mkdir -p ${output_dir}
cd ${output_dir}

if [[ ${pop} == "All" ]]; then
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_1_2male_filtered.vcf --window-pi 10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_1 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_2_2male_filtered.vcf --window-pi 10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_2 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_3_2male_filtered.vcf --window-pi 10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_3 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_4_2male_filtered.vcf --window-pi 10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_4 
else
  sqlite3 -cmd ".mode csv" -cmd '.separator "\t"' /work/users/d/t/dturissi/lutzomyia/lutzomyia_metadata.db "select sample
                                                                    from lutzomyia_samples
                                                                    where pop = '${pop}'
                                                                    union all
                                                                    select sample
                                                                    from sobral_1s_subpops
                                                                    where pop = '${pop}'
                                                                    union all
                                                                    select sample
                                                                    from jacobina_subpops
                                                                    where pop = '${pop}'" > ${pop_id_file}
  
  
  #run vcftools
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_1_2male_filtered.vcf --keep ${pop_id_file} --window-pi  10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_1 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_2_2male_filtered.vcf --keep ${pop_id_file} --window-pi  10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_2 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_3_2male_filtered.vcf --keep ${pop_id_file} --window-pi  10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_3 
  vcftools --vcf /work/users/d/t/dturissi/lutzomyia/lutzo_chr_4_2male_filtered.vcf --keep ${pop_id_file} --window-pi  10000 --out ${output_dir}/pi_vcftools_10000_${pop}_chr_4 
fi



