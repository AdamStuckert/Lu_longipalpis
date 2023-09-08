## Genome annotation

First we assembled a transcriptome with RNAseq data from Ricuarte, Cundinamarca, Colombia (same population as our population genomics). These data were provided by coauthor Dr. Rafael Vivero. Unfortunately, we only have female data. We used Trinity to assemble this transcriptome, and the [code lives here]().

We developed a pipeline to model, identify, and mask repeats. We then used these data to inform an annotation of our assemblies. Annotation was done using Maker, with our transcirptome as transcript evidence and protein evidence from Uniprot, *Drosophila melanogaster*, and *Anopheles gambiae.* Our [repeat and annotation pipeline lives here]().

We called this script with the following call for the female assembly:

```bash
sh /project/stuckert/users/Stuckert/scripts/RepeatsAnnotationtest.sh -a /project/stuckert/users/Stuckert/psychodids/Lu_long_female/Lutzomyia_longipalpis_female_1.0.fa -s Lu_longipalpis -x Female -t diptera -u /project/stuckert/users/Stuckert/peptide_databases/uniprot_sprot.fasta
```

We called this script with the following call for male assembly:

```bash
sh /project/stuckert/users/Stuckert/scripts/RepeatsAnnotationtest.sh -a /project/stuckert/users/Stuckert/psychodids/Lu_long_male/Lutzomyia_longipalpis_male_1.0.fa -s Lu_longipalpis -x Male -t diptera -u /project/stuckert/users/Stuckert/peptide_databases/uniprot_sprot.fasta
```
