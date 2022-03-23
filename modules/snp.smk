#########################
######  SNPs RULES ######
#########################



# TODO: Clair needs a lot of configuration to run not just environement



#### CLAIR #######
##################

# CLAIR Parameters
#=================

if config['read_type'] == "ccs":
    training_data=config["training_data_ccs"]
    platform="hifi"
    training_data="/bin/models/hifi"
elif config['read_type'] == "ont":
    training_data=config["training_data_ont"]
    platform="ont"
    training_data="/bin/models/ont"
elif config['read_type'] == "clr":
    platform="hifi"
    training_data=config["training_data_clr"]
    training_data="/bin/models/hifi"
else:
    print("Unknow data type, supported format are: ont, ccs, and clr")
    exit(1)

# CLAIR RULE
#===========


# CLAIR CHUNK RULE
#=================

rule callSNVsChunk:
    """
    Calling SNPs using clair
    """
    input:
        bam=data_dir + "/align/{aligner}/data.bam",
        data_index=data_dir + "/align/{aligner}/data.bam.bai",
        reference=REFERENCES,
    output:
        temp(data_dir + "/snp/{aligner}/chr.split.{chr}_{region,\d+}/merge_output.vcf.gz")
    params:
        train_data = training_data,
        platform = platform,
        start = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region)],
        end = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region) + 1]
    benchmark: data_dir + "/benchmark/snp/{aligner}/chr.split.{chr}_{region}/{chr}_{region}.benchmark.txt"
    conda: CLAIR_ENV
    log: data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/data.split.{chr}_{region}.log"
    # log: data_dir + "/snp/{aligner}/data.split.{chr,[A-Za-z0-9]+}_{region}.log"
    threads: config['clair_threads']
    shell:
        """
        echo $'{wildcards.chr}\t{params.start}\t{params.end}' > {wildcards.chr}.{params.start}.{params.end}.bed  &&\
        run_clair3.sh \
        --bam_fn {input.bam} \
        --ref_fn {input.reference} \
        --threads {threads} \
        --platform {params.platform} \
        --model_path $CONDA_PREFIX{params.train_data} \
        --output $PWD/snp/{wildcards.aligner}/chr.split.{wildcards.chr}_{wildcards.region} \
        --bed_fn={wildcards.chr}.{params.start}.{params.end}.bed > {log} 2>&1 \
        && rm {wildcards.chr}.{params.start}.{params.end}.bed
        """


#### CALL VARINAT BY CHUNKS #######
###################################


rule concatChromosome:
    """
    Concat splited chromomsomes regions
    """
    input: lambda wildcards: expand(data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/merge_output.vcf.gz", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
        #vcf_index = lambda wildcards: expand(data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/merge_output.vcf.gz.tbi", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
    output: temp(data_dir + "/snp/{aligner}/data.{chr}.vcf")
    message: "Concat variant split per Chromosome"
    conda: PRINCESS_ENV
    benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
    params:
        temp_chr=data_dir + "/snp/{aligner}/data.{chr}_filtered.vcf",
        filter=config['filter_chrs'],
        read_type=config['read_type']
    shell:"""
        bcftools concat -o {output} {input}
        """

#### UPDATE HEADER #######
##########################

rule updateHeader:
    """
    Update the phased SNPs in phased/aligner/data.vcf
    Where the PS in header defined as Integer where it should be String.
    Result: phased/aligner/data_update_header.vcf
    Will be used in mergeParentalSNPs rule later
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

rule vcfIndex:
    """
    Index VCF file.
    """
    input: data_dir + "/{sample}.vcf.gz"
    output: data_dir + "/{sample}.vcf.gz.tbi"
    message: "Indexing vcf file {input}"
    conda: PRINCESS_ENV
    shell:"""
        tabix -p vcf {input}
        """

#### MERGING PHASED VCF FILE WITH PARENTAL SNPs ########
########################################################

rule mergeParentalSNPs:
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

rule updateSNPs:
    """
    Here we shall take the input from mergeParentalSNPs but we need to unzip it first.
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
rule concactSNPs:
    """
    Rule to concat the identifed SNPs this will only be called by the user
    in case if he wanted to have only SNPs
    """
    input: lambda wildcards: expand(data_dir + "/snp/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list),
    output: data_dir + "/snp/{aligner}/data.vcf"
    message: "Concat SNP files"
    benchmark: data_dir + "/benchmark/snp/{aligner}/concat_snp.vcf"
    conda: PRINCESS_ENV
    shell:"""
        vcfcat {input} | vcfstreamsort > {output}
        """
