#!/bin/bash AlignReads.SRAdata.sh
# USAGE: sh AlignReads.SRAdata.sh -o ${WORK}/longipalpis_ms/map2male -r /nas/longleaf/home/astuck/longipalpis_ms/Lutzomyia_longipalpis_male_1.0.fa -s /nas/longleaf/home/astuck/longipalpis_ms/SRAsamples.txt -j 1

while getopts o:r:s:j: option
do
case "${option}"
in
o) OUT=${OPTARG};;
r) ref=${OPTARG};;
s) SAMPLES=${OPTARG};;
j) PERJOB=${OPTARG};;
esac
done

lutzdatdir="/work/users/a/s/astuck/sra_data"
WORK="/work/users/a/s/astuck"
gname=$(echo $ref | sed "s/.fa//")


# load modules
ml bbmap/38.87 bwa/0.7.17 gatk/4.1.9.0 picard/2.23.4
# directory to drop in temp files
mkdir -p ${OUT}
mkdir ${OUT}/fastqs
mkdir ${OUT}/bams
mkdir ${OUT}/gvcfs
mkdir ${OUT}/logs 


# index assembly
bwa index $ref
samtools faidx $ref

# create assembly dictionary
java -jar /nas/longleaf/apps/picard/2.20.0/picard-2.20.0/picard.jar CreateSequenceDictionary R=$ref O=$gname.dict


## create file for reading job ids in order to autopopulate downstream work
> alignreads_jobs.txt

## Parse out samples

#basename -a ${lutzdatdir}/*_1.fastq | sed "s/_1.fastq//g" > SRAsamples.txt

# using only "good samples"
num_samps=$(wc -l $SAMPLES | cut -f1 -d " ")
cat $SAMPLES | cut -f1 > samples.txt

printf "There are $num_samps total samples\n\n\n"

flt=$(awk "BEGIN {print $num_samps/$PERJOB}")
num_scripts=$(echo $flt | awk '
   function ceil(valor)
   {
      return (valor == int(valor)) ? valor : int(valor)+1
   }
{
   printf "%d", ceil($1)
} ')

printf "Populating $num_scripts sbatch jobs\n\n\n"

lower=1
upper=$PERJOB

for num in $(seq 1 $num_scripts)
do


cat << EOF > lutzoSNPs_$num.job
#!/bin/bash
#SBATCH --ntasks 1
#SBATCH -t 4-0:00:00
#SBATCH --job-name=lutzoSNPs
#SBATCH --output=lutzoSNPs_%j.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=general

# load modules
ml bbmap/38.87 bwa/0.7.17 gatk/4.1.9.0 picard/2.23.4

lutzdatdir="/work/users/a/s/astuck/sra_data"
OUT=$OUT
ref=$ref 

SCR="/pine/scr/a/s/astuck"

lower=$lower
upper=$upper

EOF

cat << "EOF" >> lutzoSNPs_$num.job

samples=$(sed -n "$lower","$upper"p samples.txt)

printf "Submitting jobs on: $samples\n\n"
for sample in $samples
do

# trim reads
bbduk.sh in=${lutzdatdir}/${sample}_1.fastq out=${OUT}/fastqs/${sample}_L001_1.fastq in2=${lutzdatdir}/${sample}_2.fastq out2=${OUT}/fastqs/${sample}_L001_2.fastq

# align reads
bwa mem -t 24 $ref ${OUT}/fastqs/${sample}_L001_1.fastq ${OUT}/fastqs/${sample}_L001_2.fastq | samtools view --threads 24 -Sb -o ${OUT}/bams/${sample}_L001.bam

# add read group
echo Specifying sequencing lane as lane 1
java -jar /nas/longleaf/apps/picard/2.20.0/picard-2.20.0/picard.jar AddOrReplaceReadGroups I=${OUT}/bams/${sample}_L001.bam O=${OUT}/bams/${sample}_L001.renamed.bam RGLB=lib1 RGPL=illumina RGSM=$(echo $sample) RGPU=L001 RGID=1


	ln -s ${OUT}/bams/${sample}_L001.bam ${OUT}/bams/${sample}.merged.bam



# Clean bams
java -jar /nas/longleaf/apps/picard/2.20.0/picard-2.20.0/picard.jar CleanSam I=${OUT}/bams/${sample}.merged.bam O=${OUT}/bams/${sample}.clean.bam QUIET=true COMPRESSION_LEVEL=0
# sort bams
java -jar /nas/longleaf/apps/picard/2.20.0/picard-2.20.0/picard.jar SortSam I=${OUT}/bams/${sample}.clean.bam O=${OUT}/bams/${sample}.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
# dedup BAMS
java -jar /nas/longleaf/apps/picard/2.20.0/picard-2.20.0/picard.jar MarkDuplicates I=${OUT}/bams/${sample}.sorted.bam O=${OUT}/bams/${sample}.dedupd.bam ASSUME_SORTED=true CREATE_INDEX=true M=${OUT}/logs/${sample}.dedup.metric.txt


# HaplotypeCaller: call haplotypes
gatk --java-options "-Xmx12g" HaplotypeCaller -R ${ref} -I ${OUT}/bams/${sample}.dedupd.bam -O ${OUT}/gvcfs/${sample}.gvcf -ERC GVCF

# Calculate mapping rate
mkdir ${OUT}/mapping
echo calculating mapping for ${sample}
percent=$(samtools flagstat --threads 4 ${OUT}/bams/${sample}.dedupd.bam | awk -F "[(|%]" 'NR == 7 {print $2}')
printf "$sample\t$percent\n" >> ${OUT}/mapping/mapping.tsv

## average depth (per bp)
mkdir ${OUT}/depth
echo calculating average depth for ${sample}
depth=$(samtools depth --threads 4 -a ${OUT}/bams/${sample}.dedupd.bam |  awk '{sum+=$3} END { print sum/NR}')
sd=$( samtools depth --threads 4 -a ${OUT}/bams/${sample}.dedupd.bam |  awk '{sum+=$3; sumsq+=$3*$3} END { print sqrt(sumsq/NR - (sum/NR)**2)}')
printf "$sample\t$depth\t$sd\n" >> ${OUT}/depth/seqdepth.tsv

## Per bp Depth
samtools depth --threads 4 -a  ${OUT}/bams/${sample}.dedupd.bam > ${OUT}/depth/${sample}.depth.tsv



done







EOF

# revise uppers and lowers for next iteration
lower=$(($lower + $PERJOB))
upper=$(($upper + $PERJOB))

# submit the job
sbatch lutzoSNPs_${num}.job | tee -a alignreads_jobs.txt

done



#
### Downstream analyses....
 # note.... $ sbatch --dependency=afterok:<jobID_A:jobID_C:jobID_D> jobB.sh
 # prep depdencies
cat alignreads_jobs.txt | cut -f 4  -d " " > alignmentjobs.txt

count=1
printf "#SBATCH --dependency=afterok:" > dep.txt
while read line
do
if [ $count -lt $num_scripts ]
  then
  printf "$line:" >> dep.txt
  count=$(($count + 1))
else
  printf "$line" >> dep.txt
fi
done < alignmentjobs.txt

deps=$(cat dep.txt)

############### Genomics DB Import + GenotypeGVCFs


#### change this if your samples file is a single line with location of fastqs
# prep samples file:
cat $SAMPLES | cut -f1 >  db.samples.txt
sed -i "s&^&${OUT}/gvcfs/&g" db.samples.txt

## Script header
cat << EOF > lutzoSNPS_Genotype.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=24
#SBATCH -t 5-0
#SBATCH --mem=30g
#SBATCH -o Genotype_%j.out
$deps

module purge
ml gatk/4.1.9.0

EOF

grep ">" $ref > contigs.txt
sed -i "s/>//g" contigs.txt

# make directory for db:
mkdir ${OUT}/lutzo_db
sh $HOME/longipalpis_ms/scripts/gvcf2vcf.sh -r $ref -n ${OUT}/lutzo_db -o ${OUT}/lutzo_all -c contigs.txt -s db.samples.txt >> lutzoSNPS_Genotype.job
#sbatch lutzoSNPS_Genotype.job
