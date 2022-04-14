### Analytical Workbook

This file details the analytical approaches we took for this manuscript. I will reference the scripts we ran, in the order we ran them (more or less).

## Aligning read data and calling variants

To align reads to our assembly, I used our script [AlignReads.sh](https://github.com/AdamStuckert/Lu_longipalpis/blob/main/scripts/AlignReads.sh). This aligns reads with bwa and then follows GATK best practices to clean, sort, deduplicate, and ultimately call genotypes. The script also splits samples up into a number of jobs to speed up mapping. I used this call:

>sh AlignReads.sh -o ${SCR}/All_Lutzo_male -r $PROJ/genomes/Lu_longipalpis_male_1.0.fa -s $HOME/longipalpis_ms/Samples5xCoverage.tsv -j 5
