## Analytical Workbook

This file details the analytical approaches we took for this manuscript. I will reference the scripts we ran, in the order we ran them (more or less).

### Aligning read data and calling variants

To align reads to our assembly, I used our script [AlignReads.sh](https://github.com/AdamStuckert/Lu_longipalpis/blob/main/scripts/AlignReads.sh). This aligns reads with bwa and then follows GATK best practices to clean, sort, deduplicate, and ultimately call genotypes. The script also splits samples up into a number of jobs to speed up mapping. I used this call:

>sh AlignReads.sh -o ${SCR}/All_Lutzo_male -r $PROJ/genomes/Lu_longipalpis_male_1.0.fa -s $HOME/longipalpis_ms/Samples5xCoverage.tsv -j 5

### Figure 1

This figure combines results from Principle Components Analyses (PCA) with that of Fst and Pi analyses. I created a PCA in three ways: first with the whole genome data, second with just the presumptive sex chromosome (chr 1), and finally using just autosomal data (chrs 2-4). The input data was mapped reads and I used PCAngsd for this analysis using the script [angsdPCA.job](https://github.com/AdamStuckert/Lu_longipalpis/blob/main/scripts/angsdPCA.job).

We ran Pixy to calculate both Fst and Pi. Input data was a filtered, vcf with both invariant and variant sites. I calculated Fst between males and females, and then calculated Pi for both males and females indpendently. This analysis was run with the script **ADD SCRIPT**.

Finally, the data from PCAngsd and Pixy are read into R and used to make [Figure 1 for the manuscript using this R script](https://github.com/AdamStuckert/Lu_longipalpis/blob/main/Figures/Fig1_PCA_Manhattans.Rmd).


### Phylogenetic analyses

1. Ran a full tree approach across the continent with `iqtree.sh`. 
2. Extracted BUSCO genes and got their sequences.
