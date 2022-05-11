##########################
######  ALIGN RULES ######
#########################



#### MINIMAP2 ####
##################

# Minimap2 Parameters
#====================

if config["read_type"].lower() in ["clr", "ccs"]:
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
        datain=data_dir + "/{sample}"
    output:
        dataout=temp(data_dir + "/align/minimap/{sample}.bam")
    params:
        reference=REFERENCES,
        h = minimap2_read_type,
        md = "--MD",
        x = x_param,
        # rg = "@RG\\tSM:SAMPLE\\tID:LONG", shuld be used like -R {params.rg}
    log:
        data_dir + "/align/minimap/{sample}.log"
    message:
        "Running minimap2 , sample is: {wildcards.sample}"
    threads: config['aligner_threads']
    benchmark: data_dir + "/benchmark/align/{sample}.minimap.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            minimap2  {params.x} -Y -R "@RG\\tID:hg2" "{params.reference}"  "{input.datain}" {params.h} "{params.md}" -t "{threads}" | samtools sort -@ {threads} - > "{output.dataout}"
            """

#### NGMLR ####
###############

rule ngmlr:
    """
    Using ngmlr to align reads
    """
    input:
        datain=data_dir + "/{sample}"
    output:
        dataout=temp(data_dir + "/align/ngmlr/{sample}.sam")
    params:
        reference=REFERENCES,
        platform="pacbio" if config["read_type"] in ["clr", "ccs"] else "ont" if config["read_type"] == "ont" else ""
    log:
        data_dir + "/align/ngmlr/{sample}.log"
    message:
        "Running ngmlr , sample is: {wildcards.sample}"
    threads: config['aligner_threads']
    benchmark: data_dir + "/benchmark/align/{sample}.ngmlr.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
            ngmlr -r "{params.reference}" -q "{input.datain}" --rg-sm SAMPLE -o "{output.dataout}" -t "{threads}" -x "{params.platform}" --bam-fix > {log} 2>&1
            """

#### SAM2BAM ####
################

rule sam2bam:
    input: data_dir + "/align/ngmlr/{sample}.sam"
    output: temp(data_dir + "/align/ngmlr/{sample}.bam")
    message: "Covert SAM to sorted BAM"
    threads: config['aligner_threads']
    benchmark: data_dir + "/benchmark/align/{sample}.sam2bam.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        samtools view -bhS {input} | samtools sort -@ {threads} - > {output}
        """

#### INDEX BAM ####
###################

rule indexBam:
    """
    Indexing bam file.
    """
    input:
        data_dir + "/{sample}.bam"
    output:
        temp(data_dir + "/{sample}.bam.bai")
    benchmark: data_dir + "/benchmark/align/{sample}.index.benchmark.txt"
    message: "Indexing {input}"
    conda: PRINCESS_ENV
    shell:
        "samtools index {input}"

#### MERGE BAM FILES ####
########################

rule mergeAlign:
    input:
        bams=lambda wildcards: expand(data_dir + "/align/{aligner}/{sample}.bam", aligner=wildcards.aligner, sample=sample_list),
        index_bams=lambda wildcards: expand(data_dir + "/align/{aligner}/{sample}.bam.bai", aligner=wildcards.aligner, sample=sample_list),
    output:
        file_name=temp(data_dir + "/align/{aligner}/data.bam")
    message:"Mergeing data"
    threads: config['samtools_threads']
    benchmark: data_dir + "/benchmark/align/{aligner}.merging.benchmark.txt"
    log:
        data_dir + "/align/{aligner}/merge.log"
    conda: PRINCESS_ENV
    threads: config['aligner_threads']
    shell:"""
        samtools merge -@ {threads} {output} {input.bams} > {log} 2>&1
        """

#### ADD RG TO BAM FILE ####
############################

rule addRG:
    input:data_dir + "/{sample}.bam"
    output:temp(data_dir + "/{sample}_rg.bam")
    params:
        rg = "@RG\\tSM:SAMPLE\\tID:LONG",
    conda: PRINCESS_ENV
    shell:"""
        samtools addreplacerg -r "{params.rg}" -o {output} {input}
    """



#### CONVER BAM FILE TO TAB ####
################################

rule bam2tab:
    """
    This rules takes bam file and extract to tab delimeted file: reads   HP  PS.
    """
    input:
        bam_file = data_dir + "/align/{aligner}/data_hap.bam",
    output: data_dir + "/align/{aligner}/data_hap.tab",
    message: "Extracting read hp and ps info from tagged bam file."
    conda: PRINCESS_ENV
    benchmark: data_dir + "/benchmark/align/{aligner}.bam2tab.benchmark.txt"
    shell:"""
        samtools index {input} && samtools view  {input.bam_file} |  grep  "PS:i:" |  awk 'BEGIN{{OFS="\\t";}}{{print $1,$(NF-2), $(NF)}}'  > {output}
        """
