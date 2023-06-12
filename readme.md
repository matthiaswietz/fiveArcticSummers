## Five summers across Fram Strait, Arctic Ocean: regional and vertical patterns in microbial communities and substrate regimes 

Bioinformatic code for bacterial community analyses in "Five summers across Fram Strait, Arctic Ocean: regional and vertical patterns in microbial communities and substrate regimes" by Wietz and colleagues (doi: XXX). Metabarcoding reads were obtained using 16S rRNA primers, followed by Illumina MiSeq sequencing and generation of amplicon sequence variants (ASVs) using DADA2 following https://benjjneb.github.io/dada2/tutorial_1_8.html. The repo includes the full pipeline (raw fastq reads to ASV taxonomy) and all subsequent analyses.

### Overview of content

Directory "dada":
Rmarkdowns describing primer clipping and DADA2 pipeline

Resulting [ASV table](seqtab.txt), [taxonomy table](tax.txt), and [ASV sequences](ASV_seqs.fasta)

[ENA accession numbers](ENA_accessions.txt) of all raw fastq files

[Physicochemical data and sample info](metadata.txt), needed for community analyses

Directory "analysis-plotting"
Scripts for analyses/figures, to be run in this order: [load/format data](DataLoad.R), [calculate alpha-diversity](RarefacDiversity.R), [community analyses](Res_Communities.R), [environmental correlations](Res_Correlations.R)
