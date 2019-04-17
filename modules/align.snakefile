# Minimap2 Parameters
#####################
if config["read_type"].lower() == "pacbio":
    minimap2_read_type = "-H"
    x_param = "-ax map-pb"
elif config["read_type"].lower() == "ont":
    x_param = "-ax map-ont"
    minimap2_read_type = ""
else:
    minimap2_read_type = ""
    x_param = ""
#####################


rule minimap2:
    """
    Using Minimap2 to align reads
    """
    input:
        datain = os.path.join(config['sample_directory'], "{sample}."+extension)
    output:
        dataout="align/minimap/{sample}.bam"
    params:
        reference=config['ref'], #lambda wildcards:ref[wildcards.ref]["db"],
        platform=config["read_type"],
        h = minimap2_read_type,  # either use -H in pacbio case or not
        md = "--MD",
        x = x_param,
    benchmark:
        "align/minimap/{sample}.minimap.benchmark.txt"
    conda:
        ALIGN
    message:
        "Running minimap2 , sample is: {wildcards.sample}"
    threads:config['minimap_threads']
    run:
        shell("""
            minimap2 {params.x} "{params.reference}"  "{input.datain}" {params.h} "{params.md}" -t "{threads}" | samtools sort -@ {threads} - > "{output.dataout}"
            """)

rule index_bam:
    """
    Indexing bam file.
    """
    input:
        "align/{aligner}/{sample}.bam"
    output:
        "align/{aligner}/{sample}.bam.bai"
    conda:
        ALIGN
    shell:
        "samtools index {input}"

rule merge_align:
    input:
        bams=lambda wildcards: expand("align/{aligner}/{sample}.bam", aligner=wildcards.aligner, sample=sample_list),
        index_bams=lambda wildcards: expand("align/{aligner}/{sample}.bam.bai", aligner=wildcards.aligner, sample=sample_list),
    output:
        file_name="align/{aligner}/data.bam"
    message:"mergeing data"
    conda:
        ALIGN
    threads: config['samtools_threads']
    run:
        mybam=" ".join(["-in {}".format(i) for i in input.bams])
        print("mergeing data {}".format(mybam))
        #shell("bamtools merge " + mybam + " -out {output.file_name}")
        shell(
        """
        samtools merge {output} {input}
        """
        )
