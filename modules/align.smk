##########################
######  ALIGN RULES ######
#########################



#### MINIMAP2 ####
##################

# Minimap2 Parameters
#====================

if config["read_type"].lower() == "pacbio":
    minimap2_read_type = "-H"
    x_param = "-ax map-pb"
elif config["read_type"].lower() == "ont":
    x_param = "-ax map-ont"
    minimap2_read_type = ""
else:
    minimap2_read_type = ""
    x_param = ""

rule minimap2:
    """
    Using Minimap2 to align reads
    """
    input:
        datain="{sample}"
    output:
        dataout="align/minimap/{sample}.bam"
    params:
        reference=REFERENCES[ref[0]],
        platform=config["read_type"],
        h = minimap2_read_type,
        md = "--MD",
        x = x_param,
        # rg = "@RG\\tSM:SAMPLE\\tID:LONG", shuld be used like -R {params.rg}
    log:
        "align/minimap/{sample}.log"
    message:
        "Running minimap2 , sample is: {wildcards.sample}"
    threads: config['aligner_threads']
    benchmark: "benchmark/align/minimap/{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            minimap2  {params.x} "{params.reference}"  "{input.datain}" {params.h} "{params.md}" -t "{threads}" | samtools sort -@ {threads} - > "{output.dataout}"
            """

#### NGMLR ####
###############

rule ngmlr:
    """
    Using ngmlr to align reads
    """
    input:
        datain="{sample}"
    output:
        dataout=temp("align/ngmlr/{sample}.sam")
    params:
        reference=REFERENCES[ref[0]],
        platform=config["read_type"]
    log:
        "align/ngmlr/{sample}.log"
    message:
        "Running ngmlr , sample is: {wildcards.sample}"
    threads: config['aligner_threads']
    benchmark: "benchmark/align/ngmlr/{sample}.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            ngmlr -r "{params.reference}" -q "{input.datain}" --rg-sm SAMPLE -o "{output.dataout}" -t "{threads}" -x "{params.platform}" --bam-fix > {log} 2>&1
            """

#### SAM2BAM ####
################

rule sam2bam:
    input: "{sample}.sam"
    output: "{sample}.bam"
    message: "Covert SAM to sorted BAM"
    threads: config['aligner_threads']
    benchmark: "benchmark/align/{sample}.sam2bam.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        samtools view -bhS {input} | samtools sort -@ {threads} - > {output}
        """

#### INDEX BAM ####
###################

rule index_bam:
    """
    Indexing bam file.
    """
    input:
        "{sample}.bam"
    output:
        "{sample}.bam.bai"
    benchmark: "benchmark/align/ngmlr/{sample}.benchmark.txt"
    message: "Indexing {input}"
    conda: PRINCESS_ENV
    shell:
        "samtools index {input}"

#### MERGE BAM FILES ####
########################

rule merge_align:
    input:
        bams=lambda wildcards: expand("align/{aligner}/{sample}.bam", aligner=wildcards.aligner, sample=sample_list),
        index_bams=lambda wildcards: expand("align/{aligner}/{sample}.bam.bai", aligner=wildcards.aligner, sample=sample_list),
    output:
        file_name="align/{aligner}/data.bam"
    message:"Mergeing data"
    threads: config['samtools_threads']
    benchmark: "benchmark/align/{aligner}/merging.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        samtools merge {output} {input.bams}
        """

#### ADD RG TO BAM FILE ####
############################

rule add_rg:
    input:"{sample}.bam"
    output:"{sample}_rg.bam"
    params:
        rg = "@RG\\tSM:SAMPLE\\tID:LONG",
    conda: PRINCESS_ENV
    shell:"""
        samtools addreplacerg -r "{params.rg}" -o {output} {input}
    """



#### CONVER BAM FILE TO TAB ####
################################

rule bam_2_tab:
    """
    This rules takes bam file and extract to tab delimeted file: reads   HP  PS.
    """
    input: "align/{aligner}/data_hap.bam",
    output: "align/{aligner}/data_hap.tab",
    message: "Extracting read hp and ps info from tagged bam file."
    conda: PRINCESS_ENV
    benchmark: "benchmark/align/{aligner}/bam2tab.benchmark.txt"
    shell:"""
        samtools view  {input} |  grep  "PS:i:" |  awk 'BEGIN{{OFS="\\t";}}{{print $1,$(NF-2), $(NF)}}'  > {output}
        """
