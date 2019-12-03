#########################
######  SNPs RULES ######
#########################



# TODO: Clair needs a lot of configuration to run not just environement



#### CLAIR #######
##################

# CLAIR Parameters
#=================

CLAIR=config["clair_location"]
TRAINING_DATA = config["training_data_pacbio"]
TRAINING_DATA_ONT = config["training_data_ont"]
CLAIR=config["clair_script"]
if config['read_type'] == "pacbio":
    training_data=config["training_data_pacbio"]
elif config['read_type'] == "ont":
    training_data=config["training_data_ont"]
else:
    print("Unknow data type {} supported format are: ont and pacbio".format(config['read_type']))
    exit()

# CLAIR RULE
#===========

# rule call_snps:
#     """
#     Calling SNPs using clair
#     """
#     input:
#         bam=data_dir + "/align/{aligner}/data.bam",
#         data_index=data_dir + "/align/{aligner}/data.bam.bai",
#         reference=REFERENCES[ref[0]],
#     output:
#         data_dir + "/snp/{aligner}/data.{chr}.vcf"
#     params:
#         train_data=training_data,
#         minCoverage=config['clair_coverage'],
#         clair_location=CLAIR,
#         # process_script=process_clar,
#         temp_out=data_dir + "/{chr}_vcf.tmp",
#         mypypy=config['clair_pypy'],
#     benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
#     conda: CLAIR_ENV
#     threads: config['clair_threads']
#     shell:
#         """
#         export PATH=$PWD/bin/pypy3.5-7.0.0-linux_x86_64-portable/bin:$PATH && \
#         clair.py callVarBam \
#             --chkpnt_fn {params.train_data} \
#             --bam_fn {input.bam} \
#             --ref_fn {input.reference} \
#             --minCoverage {params.minCoverage} \
#             --ctgName {wildcards.chr} \
#             --threads {threads} --call_fn {output}
#         """

# CLAIR CHUNK RULE
#=================
# TODO: add samtools to clair-env environement
rule call_snps_chunk:
    """
    Calling SNPs using clair
    """
    input:
        bam=data_dir + "/align/{aligner}/data.bam",
        data_index=data_dir + "/align/{aligner}/data.bam.bai",
        reference=REFERENCES[ref[0]],
    output:
        temp(data_dir + "/snp/{aligner}/data.split.{chr}_{region}.vcf")
    params:
        train_data = training_data,
        minCoverage = config['clair_coverage'],
        start = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region)],
        end = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region) + 1]
    benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}_{region}.benchmark.txt"
    conda: CLAIR_ENV
    threads: config['clair_threads']
    shell:
        """
        export PATH=$PWD/bin/pypy3.5-7.0.0-linux_x86_64-portable/bin:$PATH && \
        clair.py callVarBam \
            --chkpnt_fn {params.train_data} \
            --bam_fn {input.bam} \
            --ref_fn {input.reference} \
            --minCoverage {params.minCoverage} \
            --ctgName {wildcards.chr} \
            --ctgStart {params.start} \
            --ctgEnd {params.end} \
            --threads {threads} --call_fn {output}
        """

#### CALL VARINAT BY CHUNKS #######
###################################

rule concat_chromosome:
    """
    Concat splited chromomsomes regions
    """
    input: lambda wildcards: expand(data_dir + "/snp/{aligner}/data.split.{chr}_{region}.vcf", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
    output: data_dir + "/snp/{aligner}/data.{chr,[A-Za-z0-9]+}.vcf"
    message: "Concat variant split per chromomsome"
    conda: PRINCESS_ENV
    benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
    shell:"""
        vcfcat {input} | vcfstreamsort > {output}
        """

#### UPDATE HEADER #######
##########################

rule update_header:
    """
    Update the phased SNPs in phased/aligner/data.vcf
    Where the PS in header defined as Integer where it should be String.
    Result: phased/aligner/data_update_header.vcf
    Will be used in merge_parental_snps
    """
    input:data_dir + "/{sample}.vcf"
    output:data_dir + "/{sample}_update_header.vcf"
    message:"Update header file to change from float to string"
    shell:"""
        sed 's/ID=PS,Number=1,Type=Integer,Descri/ID=PS,Number=1,Type=String,Descri/' {input} > {output}
        """

#### BGZIP UPDATED FILE #######
###############################

# TODO:  The are a better way just by using {sample}.vcf as input , but debug needed cause it conflict with other rules. beside it conflict with rule bgzip_file:

# rule bgzip_vcf:
#     """
#     bgzip phased file
#     """
#     input: "{sample}/data_update_header.vcf"
#     output: "{sample}/data_update_header.vcf.gz"
#     message: "bgzip the phased file"
#     threads: config['bgzip_threads']
#     shell:"""
#         bgzip -c -@ {threads} {input} > {output}
#         """



#### INDEXING VCF FILE ########
###############################

rule vcf_index:
    """
    Index VCF file.
    """
    input: data_dir + "/{sample}.vcf.gz"
    output: data_dir + "/{sample}.vcf.gz.tbi"
    message: "Indexing phacsed vcf file"
    conda: PRINCESS_ENV
    shell:"""
        tabix -p vcf {input}
        """

#### MERGING PHASED VCF FILE WITH PARENTAL SNPs ########
########################################################

rule merge_parental_snps:
    """
    If the user wanted to update identifed SNVs this will be the first rule in sequence,
    Input: phased SNVs after updating header using update_header the bgzip and index it using
    bgzip_vcf and vcf_index respectively.
    """
    input:
        sample_snps = data_dir + "/phased/{aligner}/data_update_header.vcf.gz",
        sample_snps_index = data_dir + "/phased/{aligner}/data_update_header.vcf.gz.tbi",
        maternal_snps = config['maternal_snps'],
        paternal_snps = config['paternal_snps'],
    output: data_dir + "/phased/{aligner}/data_paternal_maternal.vcf.gz"
    message: data_dir + "/merging vcf from samplepaternal and maternal respectively"
    benchmark: data_dir + "/benchmark/snp/{aligner}/merge_parental.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        bcftools merge {input.sample_snps} {input.paternal_snps} {input.maternal_snps} | bgzip > {output}
        """

#### UPDATING PHASED SNPs ########
##################################

rule update_snps:
    """
    Here we shall take the input from merge_parental_snps but we need to unzip it first.
    """
    input: data_dir + "/phased/{aligner}/data_paternal_maternal.vcf.gz"
    output:
        updated_vcf = data_dir + "/phased/{aligner}/data_updated.vcf",
    message: "Running update SNPs"
    params:
        update_script = config['updat_snps_script'],
        phased_stat = data_dir + "/statitics/phased/phasing_stat.txt",
        block_tsv = data_dir + "/statitics/phased/blocks.tsv",
    benchmark: data_dir + "/benchmark/snp/{aligner}/update_snps.benchmark.txt"
    shell:"""
        mkdir -p statitics/phased  &&
        python {params.update_script} -i {input} -u {output.updated_vcf} -o {params.block_tsv} -s {params.phased_stat}
        """


#### CONCAT SNPs ########
#########################
rule concact_snps:
    """
    Rule to concat the identifed SNPs this will only be called by the user
    in case if he wanted to have only SNPs
    """
    input: lambda wildcards: expand(data_dir + "/snp/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list[ref[0]]),
    output: data_dir + "/snp/{aligner}/data.vcf"
    message: "Concat SNP files"
    benchmark: data_dir + "/benchmark/snp/{aligner}/concat_snp.vcf"
    conda: PRINCESS_ENV
    shell:"""
        vcfcat {input} | vcfstreamsort > {output}
        """
