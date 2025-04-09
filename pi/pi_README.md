#load samples and pops into database

python3 load_lutzomyia_samples.py

#calculate pi

sbatch get_lutzo_pi_vcftools.sbatch

#process pi results and load into database

python3 load_pi_vcftools.py

#make figures

Rscript lutzomyia_manuscript figure_4_pi.R
Rscript lutzomyia_manuscript figure_S9.R
