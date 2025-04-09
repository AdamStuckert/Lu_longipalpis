#!/bin/bash
#SBATCH -p general
#sbatch -J vcf_to_plink
#SBATCH --cpus-per-task=1
#SBATCH -t 72:00:00
#SBATCH --mem=32g
#SBATCH -o /work/users/d/t/dturissi/lutzomyia/sweep/vcf_to_plink_%j.out
#SBATCH -e /work/users/d/t/dturissi/lutzomyia/sweep/vcf_to_plink_%j.err


plink2 --vcf /proj/matutelb/data_share/longipalpis/lutzo_all_chrs_2male.vcf.gz --make-bed --maf 0.05 --geno --chr chr_1,chr_2,chr_3,chr_4 --allow-extra-chr --const-fid --out /work/users/d/t/dturissi/lutzomyia/sweep/lutzo_all_male
