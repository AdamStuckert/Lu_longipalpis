#!/bin/bash

cd longipalpis_ms/scripts/

ml fastqc
samples=$(cut -f1 ${HOME}/longipalpis_ms/Samples5xCoverage.tsv)

for sample in $samples
do
files=$(ls ${HOME}/Lutzo_data/${sample}*)
for file in $files
do
fastqc -o fastqc_results ${file}
done
done

conda activate base # contains multiqc
cd fastqc_results
multiqc .
