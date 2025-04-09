#build index 

/proj/matutelb/software/genmap-build/bin/genmap index -F /work/users/d/t/dturissi/lutzomyia/msmc2/genomes/Lutzomyia_longipalpis_male_1.0.fa -I /work/users/d/t/dturissi/lutzomyia/msmc2/genmap_index

#make mappability map

/proj/matutelb/software/genmap-build/bin/genmap map -K 150 -E 2 -I /work/users/d/t/dturissi/lutzomyia/msmc2/genmap_index -O /work/users/d/t/dturissi/lutzomyia/msmc2/genmap_map/genmap_map -t -w -bg


#prepare and index ref genome fasta

cd /work/users/d/t/dturissi/lutzomyia/msmc2/genomes
ln -s /proj/matutelb/data_share/longipalpis/genomes/Lutzomyia_longipalpis_male_1.0.fa .
bwa index Lutzomyia_longipalpis_male_1.0.fa


#make msmc2 input files

sbatch prep_msmc2.sbatch

#run mscm2

sbatch prep_msmc2.sbatch

#process msmc2 results and load into database

python3 load_msmc2_results.py

#make mansucripot plot

Rscript msmc2_plots_manuscript.R
