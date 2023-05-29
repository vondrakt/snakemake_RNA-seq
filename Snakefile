configfile: '/home/tihana/Downloads/snakemake_RNA-seq-main/config.yml'

rule all:
    input: 'run_alignment.done'

rule concatenating_files:
    input:
        pair1 = expand('{sample}', sample=config["samples_1P"]),
        pair2 = expand('{sample}', sample=config["samples_2P"])
    output:
        out1 = config['path']+config['ID']+'_R1.fastq',
        out2 = config['path']+config['ID']+'_R2.fastq'

    shell:
        'zcat {input.pair1} > {output.out1} ; zcat {input.pair2} > {output.out2}'

rule gzip_files:
    input:
        file1 = rules.concatenating_files.output.out1,
        file2 = rules.concatenating_files.output.out2
    output:
        out1 = config['path']+config['ID']+'_R1.fastq.gz',
        out2 = config['path']+config['ID']+'_R2.fastq.gz'
    shell:
        'gzip -c {input.file1} > {output.out1}; gzip -c {input.file2} > {output.out2}'

rule trim_reads:
    input:
        file1 = rules.gzip_files.output.out1,
        file2 = rules.gzip_files.output.out2,
        adapters = config['adapters']
    output:
        out = expand('{path}{ID}_{outs}.fastq.gz', path=config['path'], ID=config['ID'], outs = ['1P', '1U', '2P', '2U'])
    shell:
        'trimmomatic PE -threads 30 -phred33 {input.file1} {input.file2} {output.out} ILLUMINACLIP:{input.adapters}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'

rule quality_control:
    input:
        files = expand('{path}{ID}_{pairs}.fastq.gz', path=config['path'], ID=config['ID'], pairs=['1P', '2P'])
    shell:
        'fastqc {input.files}'

rule index_genome:
    input:
        genome_fasta = config['genome_assembly'],
        genome_annotation = config['genome_annotation']
    params:
        genome_directory = expand('{path}Genome_directory', path=config['path'])
    shell:
        'STAR --runMode genomeGenerate --runThreadN 10 --genomeDir {params.genome_directory} '
        '--sjdbGTFfile {input.genome_annotation} --genomeFastaFiles {input.genome_fasta}'

rule run_alignment:
    input:
        reads1 = config['path']+config['ID']+'_1P.fastq.gz',
        reads2 = config['path']+config['ID']+'_2P.fastq.gz'
    output: touch('run_alignment.done')
    params:
        genome_directory = rules.index_genome.params.genome_directory,
        prefix1 = config['ID']+'_1P.sorted.bam',
        prefix2 = config['ID']+'_2P.sorted.bam'
    shell:
        'STAR --genomeDir {params.genome_directory} --runThreadN 10 --readFilesIn {input.reads1} '
        '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix {params.prefix1} --outSAMtype BAM SortedByCoordinate ; '
        'STAR --genomeDir {params.genome_directory} --runThreadN 10 --readFilesIn {input.reads2} '
        '--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix {params.prefix2} --outSAMtype BAM SortedByCoordinate'
