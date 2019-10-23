#################################
######  METHYLATION RULES #######
#################################



#### NANOPOLISH INDEX ####
##########################

rule nano_index:
    """
    Preparing index to links read ids with their signal-level data in the FAST5 files
    """
    input:
        fastq_file="{sample}",
    output: "{sample}.index.readdb"
    message: "Input file is {wildcards.sample}"
    params:
        fast5_dir = lambda wildcards: ont_sample_dir[wildcards.sample]
    benchmark: "benchmark/methylation/index.{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            nanopolish index -d {params.fast5_dir} {input.fastq_file}
            """

#### NANOPOLISH METHYLATION ####
################################

rule cal_meth:
    """
    Calling Methylation
    """
    input:
        fastq_file="{sample}",
        bam_file="align/{aligner}/{sample}.bam",
        bam_index="align/{aligner}/{sample}.bam.bai",
        fastq_index="{sample}.index.readdb",
    output: "meth/{aligner}/{sample}.methylation_calls.tsv"
    params:
        ref = REFERENCES[ref[0]],
    benchmark: "benchmark/methylation/{aligner}/call_methylation.{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        nanopolish call-methylation -t 8 -r {input.fastq_file} -b {input.bam_file} -g {params.ref} > {output}
        """

#### CALL ALL METHYLATION ####
##############################

rule all_methylation:
    """
    Call all methylation samples.
    """
    input: lambda wildcards: expand("meth/{aligner}/{sample}.methylation_calls.tsv", aligner=wildcards.aligner, sample=[i for i in list(ont_sample_dir.keys())])
    output: "meth/{aligner}/methylation_calls.tsv",
    message: "Collecting all methylation samples {input}"
    shell:"""
        touch {output}
        """
