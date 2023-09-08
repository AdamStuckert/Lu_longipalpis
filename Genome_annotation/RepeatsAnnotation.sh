#!/bin/bash RepeatsAnnotation.sh
# USAGE: sh RepeatsAnnotation.sh -g GENOME -s SPECIES -x SEX -u UNIPROT
# USAGE: sh /project/stuckert/users/Stuckert/scripts/RepeatsAnnotation.sh -a /project/stuckert/users/Stuckert/R_imi_HiFi/R_imi_striped.hifiasm.bp.p_ctg.fasta -s R_imi -x Male -t vertebrata -u /project/stuckert/users/Stuckert/peptide_databases/uniprot_sprot.fasta
# USAGE: sh /project/stuckert/users/Stuckert/scripts/RepeatsAnnotationtest.sh -a /project/stuckert/users/Stuckert/psychodids/Lu_long_male/Lutzomyia_longipalpis_male_1.0.fa -s Lu_longipalpis -x Male -t diptera -u /project/stuckert/users/Stuckert/peptide_databases/uniprot_sprot.fasta

while getopts a:s:x:t:u: option
do
case "${option}"
in
a) ASSEMBLY=${OPTARG};;
s) SPECIES=${OPTARG};;
x) SEX=${OPTARG};;
t) TAXONOMY=${OPTARG};;
u) UNIPROT=${OPTARG};;
esac
done


##### NOTES: ADD IN FLAG FOR LIBRARY TAXONOMY FOR THE FUTURE
####  I NEED TO ADD IN DEPENDENCIES SO IT DOESNT SUBMIT ALL THE JOBS AT ONCE

ASSNAME=$(basename $ASSEMBLY)
SPP=$(echo $SPECIES)
PREFIX="RIMI"
#ESTS="/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_AC_3_ORP.longest_orfs.cds,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_CA_8_ORP.longest_orfs.cds,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_developmental_transcriptome.longest_orfs.cds"   #,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_PO_3_ORP.longest_orfs.cds,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_RP_10_ORP.longest_orfs.cds,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_SA_1_ORP.longest_orfs.cds,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/transcript_evidence/Ranitomeya_imitator_SR_6_ORP.longest_orfs.cds"
#PROTS="/project/stuckert/users/Stuckert/peptide_databases/uniprot_sprot.fasta,/project/stuckert/users/Stuckert/R_imi_HiFi/maker_data/protein_evidence/GCF_017654675.1_Xenopus_laevis_v10.1_protein.faa"
ESTS="/project/stuckert/users/Stuckert/psychodids/transcriptome/Trinity.fasta.transdecoder.cds.testheader.fa"
PROTS="/project/stuckert/peptide_databases/uniprot_sprot.fasta,/project/stuckert/peptide_databases/dipterans/Drosophila_melanogaster.BDGP6.32.pep.all.fa,/project/stuckert/peptide_databases/dipterans/Anopheles_gambiae.AgamP4.pep.all.fa"

# create relevant output
DIR=$(pwd)
OUT="${DIR}/${SPP}.${SEX}.annotation"
mkdir $OUT
GENOME="${OUT}/${ASSNAME}"


# symlink genome assembly
cd $OUT
ln -s $ASSEMBLY .

# tabula rasa
module purge
# load conda source script!
. /project/stuckert/software/anaconda3/etc/profile.d/conda.sh

# setup Maker variables:
UNIPROT_DB=$(echo $UNIPROT | sed "s/.fasta\$//" | sed "s/.fa\$//")



printf "Annotating our assembly $GENOME\n"
printf "This is the assembly for $SEX $SPP\n\n\n"

printf "############################################################################\n"
printf "############################################################################\n"
printf "############################################################################\n"
printf "PLEASE NOTE THAT GZIPPED GENOME FILES WILL BREAK AT THE REPEATMASKING STEP\n"
printf "############################################################################\n"
printf "I could fix this but for now I am being too lazy. \n"
printf "############################################################################\n"
printf "############################################################################\n"
printf "############################################################################\n"
# first, model repeats within the assembly
# have repeats been modeled?
if [ -f ${OUT}/${SPP}.${SEX}.repeatmodeler_db-families.fa ]
then
  printf "Repeat Modeling already completed\n\n"
else
  printf "Beginning to model repeats with RepeatModeler2\n\n"

  ## Script header
cat << EOF > ${SPP}.${SEX}.repeatmodeler.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=48
#SBATCH -t 5-0
#SBATCH --mem=100g
#SBATCH -o ${SPP}.${SEX}.repeatmodeler_%j.out
#SBATCH --mail-user=astuckert@uh.edu

module add Maker/3.01.04

printf "Running RepeatModeler2\n\n\n"


BuildDatabase -name ${SPP}.${SEX}.repeatmodeler_db -engine ncbi $GENOME
RepeatModeler -database ${SPP}.${SEX}.repeatmodeler_db -threads 48

EOF

  sbatch ${SPP}.${SEX}.repeatmodeler.job | tee repeatmodeler.sbatch.jobid.txt

  # get dependencies
  REPMODJOBID=$(tail -n1 repeatmodeler.sbatch.jobid.txt | cut -f 4 -d " ")
  rm repeatmodeler.sbatch.jobid.txt
fi

#### Taxa specific repeat library:
if [ -f ${OUT}/${TAXONOMY}.fa ]
then

  printf "$TAXONOMY repeat library exists \n\n\n"
else
  printf "Creating $TAXONOMY library\n\n\n"

ml Maker
famdb.py -i /project/dsi/apps/easybuild/software/RepeatMasker/4.1.4-foss-2020b/Libraries/RepeatMaskerLib.h5  families \
--format fasta_name --ancestors --descendants "${TAXONOMY}" --include-class-in-name > ${OUT}/${TAXONOMY}.fa

fi

#### Identify modeled repeats....
#
if [ -f ${OUT}/${TAXONOMY}.${SPP}.${SEX}.IDrepeats.fa ]
then
  printf "Modeled repeats identiied\n\n\n"
else

  printf "Identifying repeats from RepeatModeler2\n\n\n"

  ## Script header
cat << EOF > ${SPP}.${SEX}.RepeatID.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=24
#SBATCH -t 5-0
#SBATCH --mem=100g
#SBATCH -o ${SPP}.${SEX}.repeatID_%j.out
#SBATCH --mail-user=astuckert@uh.edu
#SBATCH --dependency=afterok:${REPMODJOBID}

# load conda source script!
. /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
conda activate transposon_annotation

# variables for later use in script....
SPP=$(echo $SPP)
SEX=$(echo $SEX)
OUT=$(echo $OUT)
# ID repeats:
transposon_classifier_RFSB -mode classify -fastaFile ${SPP}.${SEX}.repeatmodeler_db-families.fa -outputPredictionFile ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.txt

EOF

# need some new variables....
cat << "EOF" >> ${SPP}.${SEX}.RepeatID.job
# there are 4 lines at the end to remove!
lines=$(wc -l ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.txt | cut -f1 -d " ")
keep=$(($lines - 4))

head -n $keep ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.txt > ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.fixed.txt

# replace headers in file...
cp  ${SPP}.${SEX}.repeatmodeler_db-families.fa ${SPP}.${SEX}.repeatmodeler_IDs.fa


line=1
while [ $line -le $keep ]
do
header=$(sed -n "$line"p ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.fixed.txt)
line=$(($line + 1))
newheader=$(sed -n "$line"p ${SPP}.${SEX}.repeatmodeler_db-families.RFSB_results.fixed.txt | cut -f1 -d " ")
sed -i 's/'"${header}"'/\>'"$newheader"'/' ${SPP}.${SEX}.repeatmodeler_IDs.fa
line=$(($line + 1))
done

# now make headers into RepeatMasker format.
# now make headers into RepeatMasker format.
# format: >repeatname#class/subclass
#or simply >repeatname#class
grep "^>" ${SPP}.${SEX}.repeatmodeler_IDs.fa  | grep -v "rnd-" | sort | uniq
# this indicates there are few classes here. I'm gonna do this by hand
sed -i 's/CMC,TIR,DNATransposon/Unknown_Transib#DNA\/CMC-Transib/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Copia,LTR,Retrotransposon/Copia#LTR\/Copia/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Gypsy,LTR,Retrotransposon/Gypsy#LTR\/Gypsy/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/hAT,TIR,DNATransposon/Unknown_hAT#DNA\/hAT/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/MITE,DNATransposon/Unknown_MITE#DNA\/P/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Sola,TIR,DNATransposon/Unknown_Sola#DNA\/Sola/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Zator,TIR,DNATransposon/Unknown_Zator#DNA\/Zator/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/ERV,LTR,Retrotransposon/ERV#LTR\/ERV/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Helitron#DNA\/Helitron/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/LINE,Non-LTR/Unknown_LINE#LINE\/LINE/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Novosib,TIR,DNATransposon/Unknown_Novosib#DNA\/Novosib /g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/SINE,Non-LTR,Retrotransposon/Unknown_SINE#SINE\/SINE/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa
sed -i 's/Tc1-Mariner,TIR,DNATransposon/Tc1-Mariner#DNA\/Tc1-Mariner/g' ${SPP}.${SEX}.repeatmodeler_IDs.fa

printf "Adding %s and %s specific library to the vertebrate repeats\n\n\n" "$SPP" "$SEX"

cat ${TAXONOMY}.fa ${SPP}.${SEX}.repeatmodeler_IDs.fa  > ${TAXONOMY}.${SPP}.${SEX}.IDrepeats.fa

printf "Modeled repeats identified\n\n\n"

EOF

sbatch ${SPP}.${SEX}.RepeatID.job  | tee repeatID.sbatch.jobid.txt

  # get dependencies
  REPIDJOBID=$(tail -n1 repeatID.sbatch.jobid.txt | cut -f 4 -d " ")
  rm repeatID.sbatch.jobid.txt

fi

### Run RepeatMasker

if [ -f ${OUT}/${SPP}.${SEX}.masked.fa ]
then
  printf "Repeat masking complete\n\n\n"
else
  printf "Running Repeat Masker\n\n\n"

 ## Script header
cat << EOF > ${SPP}.${SEX}.RepeatMasker.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=48
#SBATCH -t 5-0
#SBATCH --mem=100g
#SBATCH -o ${SPP}.${SEX}.repeatMasker_%j.out
#SBATCH --mail-user=astuckert@uh.edu
#SBATCH --dependency=afterok:${REPIDJOBID}

  module add Maker/3.01.04

  RepeatMasker -pa 48 -gff -q $GENOME -lib ${OUT}/${TAXONOMY}.${SPP}.${SEX}.IDrepeats.fa

  printf "RepeatMasker done\n\n"

  printf "Preparing repeat gff3 file\n\n"

  rmOutToGFF3.pl $GENOME.out > ${SPP}.${SEX}.$SEX.prelim.gff

  cat ${SPP}.${SEX}.$SEX.prelim.gff | \
    perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \
    > ${SPP}.${SEX}.$SEX.gff

    # move masked fasta
    cp $GENOME.masked ${SPP}.${SEX}.masked.fa
EOF
sbatch ${SPP}.${SEX}.RepeatMasker.job
fi


# Prepare for genome annotation

if [ -f ${OUT}/Maker_round1a/${SPP}.${SEX}.fasta ]
then
  printf "Maker Round 1 already completed\n\n"
else
cat << EOF > ${SPP}.${SEX}.Maker1.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH -t 10-0
#SBATCH --mem=100g
#SBATCH -o ${SPP}.${SEX}.Maker1_%j.out
#SBATCH --mail-user=astuckert@uh.edu

  printf "Preparing to run Maker Round 1\n\n\n"
  mkdir Maker_round1a
  cd Maker_round1a

  # setup environment
  . /project/stuckert/software/anaconda3/etc/profile.d/conda.sh
conda deactivate
conda deactivate
  module purge
  module add Maker/3.01.04

  ## IMPORTANT: export repeatmasker libraries...
  export LIB_DIR=/project/dsi/apps/easybuild/software/RepeatMasker/4.1.4-foss-2020b/Libraries
  #export REPEATMASKER_MATRICES_DIR=/project/dsi/apps/easybuild/software/RepeatMasker/4.1.4-foss-2020b/Matrices/

  # control files for Maker
  maker -CTL
  mv maker_opts.ctl ${SPP}.${SEX}.maker_opts.ctl
 # replace est
  sed -i "s&est= #set of ESTs or assembled mRNA-seq in fasta format&est=${ESTS}&" ${SPP}.${SEX}.maker_opts.ctl ###
  sed -i "s&est2genome=0 .*$&est2genome=1&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&rmlib=.*$&rmlib=${OUT}/${TAXONOMY}.${SPP}.${SEX}.IDrepeats.fa&" ${SPP}.${SEX}.maker_opts.ctl
  # sed -i "s&rm_gff= #pre-identified repeat elements from an external GFF3 file&rm_gff=$OUT/${SPP}.${SEX}.gff&" ${SPP}.${SEX}.maker_opts.ctl
 # sed -i "s&cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)&cpus=23&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no&est2genome=1&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no&protein2genome=1&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&protein=  .*$&protein=${PROTS}&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&model_org=.*&model_org=&" ${SPP}.${SEX}.maker_opts.ctl
  sed -i "s&cpus=.*&cpus=48&" ${SPP}.${SEX}.maker_opts.ctl

  opts_ctl="${SPP}.${SEX}.maker_opts.ctl"
  bopts_ctl="maker_bopts.ctl"
  exe_ctl="maker_exe.ctl"

  echo Running maker on $genome
  echo Maker output going to $OUTPUT

  printf "Using control files: \n${opts_ctl}, \nmaker_bopts.ctl, \nmaker_exe.ctl \n\n"

maker \
  -fix_nucleotides -cpus 48 -base ${SPP}.${SEX} -quiet \
  -genome $GENOME \
  ${SPP}.${SEX}.maker_opts.ctl \
  maker_bopts.ctl \
  maker_exe.ctl

 printf "Maker annotation complete\n\n\n"
  printf "Creating fasta and gff output from maker data\n\n"

  fasta_merge -d ${SPP}.${SEX}.maker.output/${SPP}.${SEX}_master_datastore_index.log -o ${SPP}.${SEX}.fasta
EOF

# submit
sbatch ${SPP}.${SEX}.Maker1.job

fi

## create maker round 1 docs

if [ -f ${OUT}/Maker_round1a/${SPP}.${SEX}.annotations.tsv ]
then
  printf "Maker Round 1 already completed and documents made\n\n"
else
printf "Creating Maker Round 1 documents \n\n"

cat << EOF > ${SPP}.${SEX}.Maker1Docs.job
#!/bin/bash
#SBATCH -p general
#sbatch -J genotype
#SBATCH --cpus-per-task=24
#SBATCH -t 5-0
#SBATCH --mem=50g
#SBATCH -o ${SPP}.${SEX}.$SEX.Maker1_%j.out
#SBATCH --mail-user=astuckert@uh.edu

cd Maker_round1a

  # setup environment
  module purge
  module add Maker/3.01.04

  gff3_merge -d ${SPP}.${SEX}.maker.output/${SPP}_master_datastore_index.log -o ${SPP}.${SEX}.gff3 -n
  # make a lastal database
  echo running lastdb
  lastdb -p $UNIPROT_DB $UNIPROT

  echo running lastal
  lastal -P22 $UNIPROT_DB ${SPP}.${SEX}.fasta.all.maker.proteins.fasta -f BlastTab > blast.out

  echo running maker_functional_fasta
  maker_functional_fasta $UNIPROT blast.out ${SPP}.${SEX}.fasta.all.maker.proteins.fasta > ${SPP}.${SEX}.functional.proteins.fasta
  maker_functional_fasta $UNIPROT blast.out ${SPP}.${SEX}.fasta.all.maker.transcripts.fasta > ${SPP}.${SEX}.functional.transcripts.fasta
  maker_functional_gff $UNIPROT blast.out ${SPP}.${SEX}.gff3 > ${SPP}.${SEX}.functional.gff3
  maker_map_ids --prefix "$PREFIX" --justify 6 ${SPP}.${SEX}.functional.gff3 > ${SPP}.${SEX}.genome.all.id.map
  map_fasta_ids ${SPP}.${SEX}.genome.all.id.map  ${SPP}.${SEX}.functional.proteins.fasta
  map_gff_ids ${SPP}.${SEX}.genome.all.id.map  ${SPP}.${SEX}.functional.gff3
  map_fasta_ids ${SPP}.${SEX}.genome.all.id.map  ${SPP}.${SEX}.functional.transcripts.fasta

  # get annotation information for RNAseq analyses
  grep "^>" ${SPP}.${SEX}.functional.transcripts.fasta | tr -d ">" > headers.txt
  awk '{print $1}' headers.txt  > transcripts.txt
  cut -f 2 -d '"' headers.txt  | sed "s/Similar to //g" > annotations.txt
  paste transcripts.txt annotations.txt > ${SPP}.${SEX}.annotations.tsv

  printf "Maker docs created\n\n"

EOF

#sbatch ${SPP}.${SEX}.Maker1Docs.job
fi
