# Create an initial phylogenetic tree of *Leishmania* genomes in order to do full-genome alignments


We plan on doing full genome alignments of all the available *Leishmania* genomes. However, there is no good extant phylogeny for this genus. As it turns out, there are >70 genomes of various quality on NCBI. In order to make an initial estimate of a phylogeny, we downloaded all the genomes, extracted BUSCOs and built a phylogeny. We will use this as input to Cactus to do full genome alignments and create a full tree.



I downloaded genomes on 11 Nov 2022.

```bash
conda activate ncbi_datasets

# data for later:
datasets summary genome taxon leishmania --as-json-lines | dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name > leishmanig_genoes_info.tsv

#2 download genomes:
datasets download genome taxon leishmania
unzip ncbi_datasets.zip

# mv them all

mkdir genomes
mv */*/*fna genomes/
```


Next, we ran BUSCO5 on all these genomes with the Euglenozoa db10 to find conserved genes.

```bash
for file in *.fna; do

genome=$(echo $file | sed "s/.fna//g")

cat << EOF > $genome.busco.job
#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --cpus-per-task=6
#SBATCH -o $genome.busco_%j.out


# load conda source script!
. ${PROJ}/software/anaconda3/etc/profile.d/conda.sh

conda activate busco5

busco -i ${genome}.fna -f -m geno --cpu 24 -o ${genome}_out -l /proj/matutelb/users/stuckert/Leishmania_genomes/busco_downloads/lineages/euglenozoa_odb10 --offline

EOF

sbatch $genome.busco.job
done
```

Once BUSCO finished, we extracted those genes present in single copies from each assembly.


```bash

# extract amino acid seqs of complete single copy genes...
mkdir -p busco_aa
mkdir -p busco_nt
for dir in $(find . -type d -name "single_copy_busco_sequences")
do
#name=$(basename $(dirname $dir)|cut -f 1-3 -d "_" )
name=$(echo $dir |cut -f 1-3 -d "_" | cut -f 2 -d /)
echo $name
for file in ${dir}/*.faa; do
  echo $file
  filename=$(basename $file)
  echo $filename
  cp $file busco_aa/${name}_${filename}
  sed -i 's/^>/>'${name}'|/g' busco_aa/${name}_${filename} 
  cut -f 1 -d ":" busco_aa/${name}_${filename} | tr '[:lower:]' '[:upper:]' > busco_aa/${name}_${filename}.1
  mv busco_aa/${name}_${filename}.1 busco_aa/${name}_${filename}
  done
done
```


We can't build a phylogeny with genes that are only complete in one species. So we need to filter those out. Realistically, having only 2 species/samples present is also not taht helpful or informative. I found that if I use 3 samples/species I recover most of the BUSCOS (72 out of 130), but if I increase this number to 5 I get only a quarter (32 out of 130). This is due to the presence of genomes that are not at all contiguous.

```bash
# First, identify all 


for file in $(find . -name "full_table.tsv")
do
echo searching $file
grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> complete_busco_ids.txt
done


# filter out buscos in fewer than 5 genomes:
# note, this is sort of arbitrary. Finding those in 5 genomes cuts the final number of buscos in half (39 vs 72 if i use >3). total of 130 complete buscos in any genome
sort complete_busco_ids.txt | uniq -c > complete_busco_ids_with_counts.txt
awk '$NF > 3 {print $2}' complete_busco_ids_with_counts.txt > final_busco_ids.txt

# put each busco id into a single file 
while read line
do
echo $line
cat busco_aa/*_${line}.faa >> ${line}_aa.fasta
done < final_busco_ids.txt
```

Do initial alignments of each gene using MAFFT.

```bash
# align using mafft
mkdir alignments
for file in *_aa.fasta
do

cat << EOF > $file.mafft.job
#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH -o $file.mafft_%j.out

ml Core/mafft/7.490

mafft --auto --amino --thread 10 $file > alignments/$file.aln

EOF

sbatch $file.mafft.job
done
```


Make a tree for each gene using RAXML

```bash
# make gene trees
cd alignments
mkdir raxml_out
for ALN in ls *aln
do

SUF=$(basename ${ALN} |cut -f 1 -d ".")

cat << EOF > $ALN.raxml.job
#!/bin/bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 4-00:00:00
#SBATCH --mem=10g
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH -o $ALN.raxml._%j.out

ml raxml/8.2.12



raxml \
    -T 10 \
    -f a \
    -m PROTGAMMAAUTO  \
    -p 12345 \
    -x 12345 \
    -# 100 \
    -s ${ALN} \
    -n raxml_out/$SUF

EOF
sbatch $ALN.raxml.job
done
```

Get the best trees, merge 'em, build an astral consensus tree of all the genes.


```bash
# concatenate all the trees together


# run astral
module purge
. ${PROJ}/software/anaconda3/etc/profile.d/conda.sh
conda activate base
java -jar ${PROJ}/software/Astral/astral.5.7.8.jar -i whole_genome_newick.contree -o whole_genome_output_consensus.newick -t 3





```


Make a figure of the tree in R for funzies. Pop this newick tree into Cactus and do full genome MSA.









