# A snakemake pipeline for a reference based RNA-seq analysis.
The pipeline is divided into two steps. Read pre-processing and alignment Snakefiles and config files can be found in the performing_alignment.zip, after which the DE calculation is performed with the Snakefile and config from calculate_DE.zip. 
The conda enviroment with the requirement can be create with DE_analysis.yml with command:

conda env create -f environment.yml

The testing read set belongs to Desmodesmus quadricauda.
The up and down regulated genes are identified using both edgeR an DESeq2 from gene and transcript sounts.
