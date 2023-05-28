configfile: 'config.yaml'

rule concatenating_files:
    input:
        pair1 = expand('/home/tihana/snakemake_RNA-seq-main/Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L{counts}_1.SAMPLE_1000.fq.gz', counts = ['1', '2', '3']),
        pair2 = expand('/home/tihana/snakemake_RNA-seq-main/Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L{counts}_2.SAMPLE_1000.fq.gz', counts = ['1', '2', '3'])
    output:
        out1 = '/home/tihana/snakemake_RNA-seq-main/Scen_control_A_R1.fastq',
        out2 = '/home/tihana/snakemake_RNA-seq-main/Scen_control_A_R2.fastq'

    shell:
        'zcat {input.pair1} > {output.out1} ; zcat {input.pair2} > {output.out2}'

rule gzip_files:
    input:
        file1 = rules.concatenating_files.output.out1,
        file2 = rules.concatenating_files.output.out2
    shell:
        'gzip {input.file1} ; gzip {input.file2}'

rule trim_reads:
    input:
        file1 = rules.concatenating_files.output.out1 + '.gz',
        file2 = rules.concatenating_files.output.out2 + '.gz',
        adapters = '/home/tihana/snakemake_RNA-seq-main/illumina_multiplex.fa'
    output:
        out = expand('/home/tihana/snakemake_RNA-seq-main/Scen_control_A_{outs}.fastq.gz', outs = ['1P', '1U', '2P', '2U'])
    shell:
        'trimmomatic PE -threads 30 -phred33 {input.file1} {input.file2} {output.out} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule quality_control:
    input:
        files = expand('/home/tihana/snakemake_RNA-seq-main/Scen_control_A_{pairs}.fastq.gz', pairs=['1P', '2P'])
    shell:
        'fastqc {input.files}'

rule index_genome:
    input:
        genome_fasta = '/home/tihana/snakemake_RNA-seq-main/Dqua.fna',
        genome_annotation = '/home/tihana/snakemake_RNA-seq-main/Desmodesmus_quadricauda.gtf',
    params:
        genome_directory = '/home/tihana/snakemake_RNA-seq-main/DESMODESMUS_GENOME'
    shell:
        'STAR --runMode genomeGenerate --runThreadN 10 --genomeDir {params.genome_directory} '
        '--sjdbGTFfile {input.genome_annotation} --genomeFastaFiles {input.genome_fasta}'

rule run_alignment:
    input:
        reads1 = '/home/tihana/snakemake_RNA-seq-main/Scen_control_A_1P.fastq.gz',
        reads2 = '/home/tihana/snakemake_RNA-seq-main/Scen_control_A_2P.fastq.gz'
    params:
        genome_directory = rules.index_genome.params.genome_directory,
        prefix1 = 'Scen_control_A_1P.sorted.bam',
        prefix2 = 'Scen_control_A_2P.sorted.bam'
    shell:
        'STAR --genomeDir {params.genome_directory} --runThreadN 10 --readFilesIn {input.reads1} '
        '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix {params.prefix1} --outSAMtype BAM SortedByCoordinate ; '
        'STAR --genomeDir {params.genome_directory} --runThreadN 10 --readFilesIn {input.reads2} '
        '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix {params.prefix2} --outSAMtype BAM SortedByCoordinate'