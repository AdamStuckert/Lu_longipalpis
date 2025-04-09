#make genomic annotation tables

python3 /proj/matutelb/projects/gwas/scripts/make_anno_db_lutzomyia.py /proj/matutelb/data_share/longipalpis/Lu_long_male.gff3 /proj/matutelb/projects/gwas/genotype_datasets/lutzomyia_male/lutzomyia_male_anno.db

#prep gwas files

sbatch lutzomyia_male_vcf_to_plink.sbatch

mkdir /proj/matutelb/projects/gwas/genotype_datasets/lutzomyia_male
cd /proj/matutelb/projects/gwas/genotype_datasets/lutzomyia_male

sed 's/chr_//g' /proj/matutelb/projects/lutzomyia/lutzo_all_male.bim > lutzomyia_male.bim
cp /proj/matutelb/projects/lutzomyia/lutzo_all_male.bed lutzomyia_male.bed
cp /proj/matutelb/projects/lutzomyia/lutzo_all_male.fam lutzomyia_male.fam

awk -F "\t" -v OFS='\t' '{print $1, $3}' /proj/matutelb/projects/lutzomyia/Lu_longipalpis_metadata.tsv > /proj/matutelb/projects/gwas/scratch/lutzomyia_sex/lutzomyia_meta_sex.txt

sed 's/Female/1/g' lutzomyia_meta_sex.txt | sed 's/Male/2/g' > lutzomyia_meta_sex_coded.txt


#run gwas

python3 /proj/matutelb/projects/gwas/scripts/run_gwas.py	/proj/matutelb/projects/gwas/scratch/lutzomyia_sex/gwas_lutzomyia_sex_male.json

Rscript /proj/matutelb/projects/gwas/scripts/gwas_snp_reports_lutzomyia_male.R lutzomyia_sex_male 1e-7


