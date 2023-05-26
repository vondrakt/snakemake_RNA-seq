rule concatenating_files:
    input:
        pair_1 = './Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L1_1.SAMPLE_1000.fq.gz ./Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L2_1.SAMPLE_1000.fq.gz ./Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L3_1.SAMPLE_1000.fq.gz',
        pair_2 = './Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L1_2.SAMPLE_1000.fq.gz ./Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L2_2.SAMPLE_1000.fq.gz ./Scen1/Scen1_EKDL220016234-1A_HN3NKDSX3_L3_2.SAMPLE_1000.fq.gz'
    output:
        out1 = 'Scen_control_A_R1.fastq',
        out2 = 'Scen_control_A_R2.fastq'
    shell:
        """
        zcat {input.pair_1} > {output.out1}
        zcat {input.pair_2} > {output.out2}
        """

# rule gzip_files:
#     input:
#         file1 = rules.concatenating_files.output.out1,
#         file2 = rules.concatenating_files.output.out2
#     output:
#         out1 = 'Scen_control_A_R1.fastq.gz',
#         out2 = 'Scen_control_A_R2.fastq.gz'
#     shell:
#         """
#         gzip {input.file1}
#         gzip {input.file2}
#         """
