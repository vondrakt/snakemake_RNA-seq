#!/usr/bin/env python3
from optparse import OptionParser
import os
import re

# THIS IS A MASTER PYTHON SCRIPT FOR RUNNING A REFERENCE BASED RNA-SEQ ANALYSIS

# defining the arguments that need to be passed to the scripts
arguments = OptionParser()

arguments.add_option('-s', '--samples', dest='samples', help='sample IDs, separated by comma. THE NAME STRUCTURE OF SAMPLES, REPLICATES AND PAIRS LOOKS LIKE SAMPLES_REPLICATES_PAIRS.FASTQ.GZ')
arguments.add_option('-r', '--replicates', dest='replicates', help='replicate IDs, separated by comma')
arguments.add_option('-p', '--pairs', dest='pairs', help='pair IDs, separated by comma')
arguments.add_option('-g', '--genome', dest='genome', help='reference genome to be used for the analysis')
arguments.add_option('-a', '--annotation-gtf', dest='annotation_gtf', help='annotation of the reference genome in gtf format')
arguments.add_option('-b', '--annotation-gff', dest='annotation_gff', help='annotation of the reference genome in gff format')
arguments.add_option('-c', '--configs', dest='config_DE', help='list of configuration files containing the binary comparisons separated by a comma. The names of config files should be in format X_vs_Y_config_DE.txt')

(options, args) = arguments.parse_args()
if options.samples is None or options.replicates is None or options.pairs is None\
        or options.genome is None or options.annotation_gtf is None or \
        options.annotation_gff is None or options.config_DE is None:
    # if one of the arguments is missing
    print('\n----------> A mandatory option is missing !\n')  # raise an error
    arguments.print_help()  # and provide the help
    exit(-1)  # exit the script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Due to incompatibility with rnaspades, this step is skipped

# Preparing the reference genome for the STAR alignment outside of the loop (so that it is not prepared separately for each replicate)
STAR_command_reference = 'STAR --runMode genomeGenerate --runThreadN 30 --genomeDir Genome_directory --sjdbGTFfile %s --genomeFastaFiles %s' % (options.annotation_gtf, options.genome)
print(STAR_command_reference)
#os.system(STAR_command_reference)

# Processing names
samples = options.samples.split(',')
replicates = options.replicates.split(',')
pairs = options.pairs.split(',')
print(samples)
print(replicates)
print(pairs)

for sample in samples:
    for replicate in replicates:
        name1 = sample+'_'+replicate+'_'+pairs[0]
        name2 = sample+'_'+replicate+'_'+pairs[1]
        print(name1, name2)

        # Trimming the reads
        trimming_command = "trimmomatic PE -threads 30 -phred33 " \
                           "%s.fastq.gz %s.fastq.gz %s_1P.fastq.gz %s_1U.fastq.gz %s_2P.fastq.gz %s_2U.fastq.gz " \
                           "ILLUMINACLIP:illumina_multiplex.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36" % (name1, name2, name1, name1, name2, name2)
        print(trimming_command)
        #os.system(trimming_command)

        # Quality checking
        fastqc_command = "fastqc %s_1P.fastq.gz %s_2P.fastq.gz" % (name1, name2)
        multiqc_command = "multiqc %s_1P.fastq.gz %s_2P.fastq.gz" % (name1, name2)
        print(fastqc_command)
        print(multiqc_command)
        #os.system(fastqc_command)
        #os.system(multiqc_command)

        # STAR alignment
        STAR_command_name1 = "STAR --genomeDir Genome_directory --runThreadN 10 --readFilesIn %s_1P.fastq.gz " \
                             "--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s_1P.sorted.bam " \
                             "--outSAMtype BAM SortedByCoordinate" % (name1, name1)
        STAR_command_name2 = "STAR --genomeDir Genome_directory --runThreadN 10 --readFilesIn %s_2P.fastq.gz " \
                             "--readFilesCommand zcat --quantMode GeneCounts --outFileNamePrefix %s_2P.sorted.bam " \
                             "--outSAMtype BAM SortedByCoordinate" % (name2, name2)
        print(STAR_command_name1)
        print(STAR_command_name2)
        #os.system(STAR_command_name1)
        #os.system(STAR_command_name2)

        # Stringtie
        stringtie_command_name1 = "stringtie %s_1P.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G %s -e -o %s.gtf -A %s.gene_abundances.tsv" % (name1, options.annotation_gff, name1, name1)
        stringtie_command_name2 = "stringtie %s_2P.sorted.bamAligned.sortedByCoord.out.bam -p 30 -G %s -e -o %s.gtf -A %s.gene_abundances.tsv" % (name2, options.annotation_gff, name2, name2)
        print(stringtie_command_name1)
        print(stringtie_command_name2)
        #os.system(stringtie_command_name1)
        #os.system(stringtie_command_name2)


# Create the samples.txt file
#print(options.config_DE)
configs = options.config_DE.split(',')
#print(configs)
for config in configs:
    #print(config)
    prepDE_name = '%s_samples.txt' % re.sub('_config_DE.txt', '', config)
    out = open(prepDE_name, 'w')
    with open(config) as c:
        for line in c:
            #print(line)
            items = line.split()
            #print(items)
            out.write(items[1]+'\t'+items[1]+'.gtf'+'\n')
    out.close()
    config_name = re.sub('_config_DE.txt', '', config)
    command = 'prepDE.py -i %s -g %s_gene_count_matrix.csv -t %s_transcript_count_matrix.csv' % (prepDE_name, config_name, config_name)
    #os.system(command)

    gene_csv = '%s_gene_count_matrix.csv' % config_name
    gene_tsv = '%s_gene_count_matrix.tsv' % config_name
    out = open(gene_tsv, 'w')
    with open(gene_csv) as g:
        for line in g:
            line = re.sub(',', '\t', line)
            out.write(line)
    out.close()

    edgeR_analysis = "run_DE_analysis.pl --matrix %s --samples_file %s " \
                     "--reference_sample condition1 --method edgeR --output %s.edgeR_gene_counts" % (
                     gene_tsv, config, config_name)
    DESeq2_DE_analysis = "run_DE_analysis.pl --matrix %s --samples_file %s " \
                         "--reference_sample condition1 --method DESeq2 --output %s.DESeq2_gene_counts" % (
                         gene_tsv, config, config_name)
    print(edgeR_analysis)
    print(DESeq2_DE_analysis)
    # os.system(edgeR_analysis)
    # os.system(DESeq2_DE_analysis)


    transcript_csv = '%s_transcript_count_matrix.csv' % config_name
    transcript_tsv = '%s_transcript_count_matrix.tsv' % config_name
    out = open(transcript_tsv, 'w')
    with open(transcript_csv) as g:
        for line in g:
            line = re.sub(',', '\t', line)
            out.write(line)
    out.close()
    edgeR_analysis = "run_DE_analysis.pl --matrix %s --samples_file %s " \
                     "--reference_sample condition1 --method edgeR --output %s.edgeR_transcript_counts" % (
        transcript_tsv, config, config_name)
    DESeq2_DE_analysis = "run_DE_analysis.pl --matrix %s --samples_file %s " \
                         "--reference_sample condition1 --method DESeq2 --output %s.DESeq2_transcript_counts" % (
        transcript_tsv, config, config_name)

    print(edgeR_analysis)
    print(DESeq2_DE_analysis)
    # os.system(edgeR_analysis)
    # os.system(DESeq2_DE_analysis)








#
#     # Extracting differentially expressed transcripts and generating heatmaps
#     # Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed
#     # at a significance of <= 0.05 in any of the pairwise sample comparisons:
#
#     # The working directory needs to be changed
#     new_dir = './edgeR_genes'
#     print(new_dir)
#     os.chdir(new_dir)
#     analyze_DE = '$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../gene_count_matrix.tsv --samples %s -P 0.05 -C 2' % options.config_DE
#     print(analyze_DE)
#     os.system(analyze_DE)
#     os.chdir('../')
#
#     new_dir = './DESeq2_genes'
#     print(new_dir)
#     os.chdir(new_dir)
#     analyze_DE = '$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../gene_count_matrix.tsv --samples %s -P 0.05 -C 2' % options.config_DE
#     print(analyze_DE)
#     os.system(analyze_DE)
#     os.chdir('../')
#
# else:
#     out_tsv = open('./transcript_count_matrix.tsv', 'w')
#     with open('./transcript_count_matrix.csv') as c:
#         for line in c:
#             line = re.sub(',', '\t', line)
#             out_tsv.write(line)
#     out_tsv.close()
#
#     DE_analysis_command_edgeR = '$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix transcript_count_matrix.tsv --samples_file %s ' \
#                                 '--reference_sample condition1 --method edgeR --output edgeR_genes' % options.config_DE
#     print(DE_analysis_command_edgeR)
#     os.system(DE_analysis_command_edgeR)
#
#     DE_analysis_command_DESeq2 = '$TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix transcript_count_matrix.tsv --samples_file %s ' \
#                                  '--reference_sample condition1 --method DESeq2 --output DESeq2_genes' % options.config_DE
#     print(DE_analysis_command_DESeq2)
#     os.system(DE_analysis_command_DESeq2)
#
#     # Extracting differentially expressed transcripts and generating heatmaps
#     # Extract those differentially expressed (DE) transcripts that are at least 4-fold differentially expressed
#     # at a significance of <= 0.05 in any of the pairwise sample comparisons:
#
#     # The working directory needs to be changed
#     new_dir = './edgeR_genes'
#     print(new_dir)
#     os.chdir(new_dir)
#     analyze_DE = '$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../transcript_count_matrix.tsv --samples %s -P 0.05 -C 2' % options.config_DE
#     print(analyze_DE)
#     os.system(analyze_DE)
#     os.chdir('../')
#
#     new_dir = './DESeq2_genes'
#     print(new_dir)
#     os.chdir(new_dir)
#     analyze_DE = '$TRINITY_HOME/Analysis/DifferentialExpression/analyze_diff_expr.pl --matrix ../transcript_count_matrix.tsv --samples %s -P 0.05 -C 2' % options.config_DE
#     print(analyze_DE)
#     os.system(analyze_DE)
#     os.chdir('../')
