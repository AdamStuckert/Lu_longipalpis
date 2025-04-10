#!/bin/bash  iqtree.sh
# USAGE: sh iqtree.sh -v /work/users/a/s/astuck/longipalpis_ms/map2male/pixy_allsites_filtered_biallelic/lutzo_all_sites.filtered.vcf.gz -g /nas/longleaf/home/astuck/longipalpis_ms/Lutzomyia_longipalpis_male_1.0.fa
# USAGE: sh iqtree.sh -v /work/users/a/s/astuck/longipalpis_ms/map2male/pixy_allsites_filtered_biallelic/lutzo_all_sites.filtered.biallelic.continent.vcf.gz -g /nas/longleaf/home/astuck/longipalpis_ms/Lutzomyia_longipalpis_male_1.0.fa

while getopts v:g: option
do
case "${option}"
in
v) VCF=${OPTARG};;
g) GENOME=${OPTARG};;
esac
done

## $1 == vcf with called/filtered genotypes
#VCF="${SCR}//longipalpis_ms/pixy_allsites_filtered/lutzo_all_sites.filtered.vcf.gz"

mkdir iqtree 
cd iqtree 

# extract lengths from chromosomes
cat $GENOME | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > $GENOME.lengths

#chr=$(awk '{print $1}' $GENOME.lengths | head -n4)
#lengths=$(awk '{print $2}' $GENOME.lengths | head -n4) # for some reason the operand is not working this way, so for now use just the 
chr=(chr_1 chr_2 chr_3 chr_4)
lengths=(49086602 36859591 29392114 22566980)

wnd=100000
for j in {0..4}; do
        c=${chr[${j}]}
        START=0
        END=${lengths[${j}]}/$wnd+1
        MAX=${lengths[${j}]}
        for ((i=$START;i<=$END;i++)); do
                        let start=$i*$wnd+1
                        let end=$start+$wnd
                        end=$((end<MAX ? end : MAX))
          
cat << "EOF" > Tree_${c}_${i}.job
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

#ml samtools iqtree

cd ${WORK}/longipalpis_ms/map2male
mkdir iqtree_continent 
cd iqtree_continent

EOF

cat << EOF >> Tree_${c}_${i}.job
##${c}_${i}
##${c}_${i}

bcftools view -r ${c}:${start}-${end} ${VCF} > ${c}_${i}_wnd${wnd}.vcf
python ~/scripts/vcf2phylip.py -i ${c}_${i}_wnd${wnd}.vcf
rm ${c}_${i}_wnd${wnd}.vcf
iqtree2 -s ${c}_${i}_wnd${wnd}.min4.phy -redo -pre ${c}_${i} -st DNA -nt 4 -ntmax 4 -bb 1000 -m MFP+ASC 

if [ -f ${c}_${i}.varsites.phy ]
then
        iqtree2 -s ${c}_${i}.varsites.phy -redo -pre ${c}_${i} -nt 4 -ntmax 4 -bb 1000 -st DNA -m MFP+ASC
else
        echo snicky snack
fi
EOF

sbatch Tree_${c}_${i}.job
done
done


###note make sure to search for failed jobs due to variant sites...
