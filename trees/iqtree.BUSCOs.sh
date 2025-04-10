#!/bin/bash

BUSCOS=/nas/longleaf/home/astuck/longipalpis_ms/scripts/full_table_busco_format.tsv
VCF=/work/users/a/s/astuck/longipalpis_ms/map2male/lutzo_all_withSRA_2male.vcf.gz

ml samtools iqtree

DIR=$(pwd)
# remove header
grep -v "Busco" $BUSCOS > tmp
# get only Complete BUSCO genes
awk 'BEGIN{FS="\t";OFS="\t"} $2== "Complete" {print $1,$3,$4,$5,$6,$10}' tmp > ${DIR}/CompleteBUSCOs.tsv
rm tmp

mkdir BUSCO_genes
cd BUSCO_genes

# extract from vcfs...
while read line
do 
ID=$(echo $line | awk '{print $1}')
GENE=$(echo $line | awk '{print $6}')

echo ID is $ID
echo gene is $GENE

# make your bed
echo $line | awk 'BEGIN { OFS = "\t"} {print $2,$3,$4}' > ${ID}.bed
CHR=$(echo $line | awk '{print $2}')
START=$(echo $line | awk '{print $3}')
END=$(echo $line | awk '{print $4}')

# write script...
cat << "EOF" > Tree_${ID}.job
#!/bin/bash
#SBATCH -p general
#sbatch -J Tree
#SBATCH --cpus-per-task=4
#SBATCH -t 5:00:00 
#SBATCH --mem=8g
#SBATCH -o trees_%j.out

module purge
#. ${PROJ}/software/anaconda3/etc/profile.d/conda.sh
#conda activate base 
ml samtools iqtree

cd ${WORK}/longipalpis_ms/map2male
mkdir iqtree_continent_BUSCOs 
cd iqtree_continent_BUSCOs

EOF

cat << EOF >> Tree_${ID}.job


bcftools view -r ${CHR}:${START}-${END} ${VCF} > ${ID}.vcf
python ~/scripts/vcf2phylip.py -i ${ID}.vcf
iqtree2 -s ${ID}.min4.phy -redo -pre ${ID}  -nt 4 -ntmax 4 -bb 1000 -st DNA 


if [ -f ${ID}.varsites.phy ]
then
        echo 
        iqtree2 -s ${ID}.varsites.phy -redo -pre ${ID} -nt 4 -ntmax 4 -bb 1000 -st DNA 
else
        echo snicky snack tree finished
fi
EOF

sbatch Tree_${ID}.job


done <  ${DIR}/CompleteBUSCOs.tsv


