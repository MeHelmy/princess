#################################
######  METHYLATION RULES #######
#################################



#### NANOPOLISH INDEX ####
##########################

rule nanoIndex:
    """
    Preparing index to links read ids with their signal-level data in the FAST5 files
    """
    input:
        fastq_file=data_dir + "/{sample}",
    output: data_dir + "/{sample}.index.readdb"
    message: "Input file is {wildcards.sample}"
    params:
        fast5_dir = config['fast5_dir']#lambda wildcards: ont_sample_dir[wildcards.sample]
    benchmark: data_dir + "/benchmark/methylation/index.{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            nanopolish index -d {params.fast5_dir} {input.fastq_file}
            """

#### NANOPOLISH METHYLATION ####
################################

rule callMeth:
    """
    Calling Methylation
    """
    input:
        fastq_file=data_dir + "/{sample}",
        bam_file=data_dir + "/align/{aligner}/{sample}.bam",
        bam_index=data_dir + "/align/{aligner}/{sample}.bam.bai",
        fastq_index=data_dir + "/{sample}.index.readdb",
    output: data_dir + "/meth/{aligner}/{sample}.methylation_calls.tsv"
    params:
        ref = REFERENCES,
    threads: config['methylation_threads']
    benchmark: data_dir + "/benchmark/methylation/{aligner}/call_methylation.{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        nanopolish call-methylation -t 8 -r {input.fastq_file} -b {input.bam_file} -g {params.ref} > {output}
        """

#### CALL ALL METHYLATION ####
##############################

rule allMethylation:
    """
    Call all methylation samples.
    """
    input: lambda wildcards: expand(data_dir + "/meth/{aligner}/{sample}.methylation_calls.tsv", aligner=wildcards.aligner, sample=config['sample_list'].split())
    output: data_dir + "/meth/{aligner}/methylation_calls.tsv",
    message: "Collecting all methylation samples {input}"
    shell:"""
        touch {output}
        """
