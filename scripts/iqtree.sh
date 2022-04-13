#!/bin/bash  iqtree.sh
# USAGE: sh iqtree.sh -v ${SCR}/longipalpis_ms/pixy_allsites_filtered/lutzo_all_sites.filtered.vcf.gz
while getopts v: option
do
case "${option}"
in
v) VCF=${OPTARG};;
esac
done

## $1 == vcf with called/filtered genotypes
VCF="${SCR}//longipalpis_ms/pixy_allsites_filtered/lutzo_all_sites.filtered.vcf.gz"

mkdir iqtree 
cd iqtree 

chr=(HiC_scaffold_1 HiC_scaffold_2 HiC_scaffold_3 HiC_scaffold_4)
lengths=(49693208 37555792 29959561 23723792)
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
#SBATCH --mem=4g
#SBATCH -o trees_%j.out

module purge
. ${PROJ}/software/anaconda3/etc/profile.d/conda.sh
conda activate base 
ml samtools/1.15 iqtree/1.6.12

cd ${SCR}/longipalpis_ms/
mkdir iqtree 
cd iqtree

EOF

cat << EOF >> Tree_${c}_${i}.job
##${c}_${i}
##${c}_${i}

bcftools view -r ${c}:${start}-${end} ${VCF} > ${c}_${i}_wnd${wnd}.vcf
python ~/scripts/vcf2phylip.py -i ${c}_${i}_wnd${wnd}.vcf
rm ${c}_${i}_wnd${wnd}.vcf
iqtree -redo -pre ${c}_${i} -st DNA -nt 4 -ntmax 4 -bb 1000 -s ${c}_${i}_wnd${wnd}.min4.phy

EOF

sbatch Tree_${c}_${i}.job
done
done
