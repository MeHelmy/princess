#########################
######  SNPs RULES ######
#########################


#### CLAIR #######
##################

# CLAIR Parameters
#=================

# if config["clair_model"]:
#     training_data=config["clair_model"]
def platform(wildcards):
    if config['read_type'] == "ccs":
        return "hifi"
    elif config['read_type'] == "ont":
        return "ont"
    elif config['read_type'] == "clr":
        return "hifi"
    else:
        print("Unknow data type, supported format are: ont, ccs, and clr")
        exit(1)

def get_model(conda_dir):
    training_data = ""
    if config['read_type'] == "ccs":
        training_data=config["clair_model"] if config["clair_model"] else None
    elif config['read_type'] == "ont":
        training_data=config["clair_model"] if config["clair_model"] else None
    elif config['read_type'] == "clr":
        training_data=config["clair_model"] if config["clair_model"] else None
    else:
        print("Unknown data type, supported format are: ont, ccs, and clr")
        exit(1)
    return training_data
# if config['read_type'] == "ccs":
#     # training_data=config["training_data_ccs"]
#     platform="hifi"
#     training_data=config["clair_model"] if config["clair_model"] else os.path.join(os.environ['CONDA_PREFIX'], "bin/models/hifi")
# elif config['read_type'] == "ont":
#     # training_data=config["training_data_ont"]
#     platform="ont"
#     training_data=config["clair_model"] if config["clair_model"] else os.path.join(os.environ['CONDA_PREFIX'], "bin/models/ont")
#     # training_data="/bin/models/ont"
# elif config['read_type'] == "clr":
#     platform="hifi"
#     # training_data=config["training_data_clr"]
#     training_data=config["clair_model"] if config["clair_model"] else os.path.join(os.environ['CONDA_PREFIX'], "bin/models/hifi")
#     # training_data="/bin/models/hifi"
# else:
#     print("Unknow data type, supported format are: ont, ccs, and clr")
#     exit(1)


# CLAIR RULE
#===========


# CLAIR CHUNK RULE
#=================

#  [ ! -f {output.gvcf} ] && cp {output.vcf} {output.gvcf} && cp {output.vcf}.tbi {output.gvcf}.tbi &&\
if config['gvcf_snv']:
    rule callSNVsChunk:
        """
        Calling SNPs using clair in case gVCF file is required.
        """
        input:
            bam=data_dir + "/align/{aligner}/data.bam",
            data_index=data_dir + "/align/{aligner}/data.bam.bai",
            reference=REFERENCES,
        output:
            vcf = temp(data_dir + "/snp/{aligner}/chr.split.{chr}_{region,\d+}/merge_output.vcf.gz"),
            gvcf = temp(data_dir + "/snp/{aligner}/chr.split.{chr}_{region,\d+}/merge_output.gvcf.gz")
        params:
            train_data = lambda wildcards: get_model(os.environ['CONDA_PREFIX']),
            platform = platform,
            start = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region)],
            end = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region) + 1],
            gvcf = "--gvcf",
        benchmark: data_dir + "/benchmark/snp/{aligner}/chr.split.{chr}_{region}/{chr}_{region}.benchmark.txt"
        conda: CLAIR_ENV
        log: data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/data.split.{chr}_{region}.log"
        threads: config['clair_threads']
        shell:
            """
            if [ {params.train_data} == "None" ]
                then
                model="$CONDA_PREFIX/bin/models/{params.platform}"
            else
                model="{params.train_data}"
            fi

            echo $'{wildcards.chr}\t{params.start}\t{params.end}' > {wildcards.chr}.{params.start}.{params.end}.bed  &&\
            run_clair3.sh \
            --bam_fn {input.bam} \
            --ref_fn {input.reference} \
            --threads {threads} \
            --platform {params.platform} \
            --model_path $model \
            --output $PWD/snp/{wildcards.aligner}/chr.split.{wildcards.chr}_{wildcards.region} \
            --bed_fn={wildcards.chr}.{params.start}.{params.end}.bed \
            {params.gvcf} > {log} 2>&1 \
            &&\
            if [ ! -f {output.gvcf} ]; then
                cp {output.vcf} {output.gvcf} 
            fi &&\
            rm {wildcards.chr}.{params.start}.{params.end}.bed
            """
else:
    rule callSNVsChunk:
        """
        Calling SNPs using clair
        """
        input:
            bam=data_dir + "/align/{aligner}/data.bam",
            data_index=data_dir + "/align/{aligner}/data.bam.bai",
            reference=REFERENCES,
        output:
            vcf = temp(data_dir + "/snp/{aligner}/chr.split.{chr}_{region,\d+}/merge_output.vcf.gz")
        params:
            train_data = lambda wildcards: get_model(os.environ['CONDA_PREFIX']),
            platform = platform,
            start = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region)],
            end = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region) + 1],
        benchmark: data_dir + "/benchmark/snp/{aligner}/chr.split.{chr}_{region}/{chr}_{region}.benchmark.txt"
        conda: CLAIR_ENV
        log: data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/data.split.{chr}_{region}.log"
        threads: config['clair_threads']
        shell:
            """
            if [ {params.train_data} == "None" ]
                then
                model="$CONDA_PREFIX/bin/models/{params.platform}"
            else
                model="{params.train_data}"
            fi

            echo $'{wildcards.chr}\t{params.start}\t{params.end}' > {wildcards.chr}.{params.start}.{params.end}.bed  &&\
            run_clair3.sh \
            --bam_fn {input.bam} \
            --ref_fn {input.reference} \
            --threads {threads} \
            --platform {params.platform} \
            --model_path $model \
            --output $PWD/snp/{wildcards.aligner}/chr.split.{wildcards.chr}_{wildcards.region} \
            --bed_fn={wildcards.chr}.{params.start}.{params.end}.bed  > {log} 2>&1 \
            && rm {wildcards.chr}.{params.start}.{params.end}.bed
            """
# resources:
    # mem_mb=lambda wildcards, attempt: 1024 * (attempt + 1) if attempt < 3
#         --model_path $CONDA_PREFIX{params.train_data} \
    # --model_path {params.train_data} \

#### CALL VARIANT BY CHUNKS #######
###################################

if config['gvcf_snv']:
    ## TODO: This function will raise an error if there is missing gvcf, missing gvcf results from calling varaint in alt contigs where there are no variants this no gVCF, solution is to check if the vcf file is empty then just copy it with gVCF name.
    # if $(head -n 1000 {input} | grep -q -v "#") ; then
    #     vcfcat {input}  | vcfstreamsort > {params.temp_chr}\
    #     && first_max=$(find_max {params.temp_chr} {params.read_type})\
    #     && threshold=$(filsn {params.temp_chr} $first_max)\
    #     && awk -v threshold=$threshold '/^#/{{print}} !/^#/{{if ( $6 >= threshold ) {{print $0}}}}' {params.temp_chr} | awk '/^#/ {{ print }} !/^#/ {{ if ($4 != $5 ) {{ print }} }}' > {output}
    # else
    #     cp {input} {output}
    # fi
    rule concatChromosome:
        """
        Concat split chromosomes regions in case of gVCF file is required.
        """
        input:
            vcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/merge_output.vcf.gz", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
            gvcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/merge_output.gvcf.gz", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
        output:
            vcf = temp(data_dir + "/snp/{aligner}/data.{chr}.vcf"),
            gvcf = temp(data_dir + "/snp/{aligner}/data.{chr}.gvcf")
        message: "Concat variant split per Chromosome"
        params:
            tmp_dir=config["tmp_directory"],
        conda: VARIANT_ENV
        benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
        shell:"""
        if [[ ! -z "{params.tmp_dir}" ]]; then
            bcftools concat {input.vcf} | bcftools sort -T {params.tmp_dir} > {output.vcf} &&\
            bcftools concat {input.gvcf} | bcftools sort -T {params.tmp_dir} > {output.gvcf}
        else
            bcftools concat {input.vcf} | bcftools sort  > {output.vcf} &&\
            bcftools concat {input.gvcf} | bcftools sort > {output.gvcf}
        fi
            """
else:
    rule concatChromosome:
        """
        Concat split chromosomes regions
        """
        input:
            vcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/chr.split.{chr}_{region}/merge_output.vcf.gz", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
        output:
            vcf = temp(data_dir + "/snp/{aligner}/data.{chr}.vcf")
        message: "Concat variant split per Chromosome"
        params:
            tmp_dir=config["tmp_directory"],
        conda: VARIANT_ENV
        benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
        shell:"""
         if [[ ! -z "{params.tmp_dir}" ]]; then
            bcftools concat {input.vcf} | bcftools sort  -T {params.tmp_dir} > {output.vcf} 
         else
             bcftools concat {input.vcf} | bcftools sort  > {output.vcf} 
         fi
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


#### INDEXING VCF FILE ########
###############################

rule vcfIndex:
    """
    Index VCF file.
    """
    input: data_dir + "/{sample}.vcf.gz"
    output: data_dir + "/{sample}.vcf.gz.tbi"
    message: "Indexing vcf file {input}"
    conda: VARIANT_ENV
    shell:"""
        tabix -p vcf {input}
        """

rule gvcfIndex:
    """
    Index VCF file.
    """
    input: data_dir + "/{sample}.gvcf.gz"
    output: data_dir + "/{sample}.gvcf.gz.tbi"
    message: "Indexing vcf file {input}"
    conda: VARIANT_ENV
    shell:"""
        tabix -p vcf {input}
        """


#### MERGING PHASED VCF FILE WITH PARENTAL SNPs ########
########################################################

rule mergeParentalSNPs:
    """
    If the user wanted to update identified SNVs this will be the first rule in sequence,
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
    conda: VARIANT_ENV
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
        phased_stat = data_dir + "/statistics/phased/phasing_stat.txt",
        block_tsv = data_dir + "/statistics/phased/blocks.tsv",
    benchmark: data_dir + "/benchmark/snp/{aligner}/update_snps.benchmark.txt"
    conda: READ_STAT_ENV
    shell:"""
        mkdir -p statistics/phased  &&
        python {params.update_script} -i {input} -u {output.updated_vcf} -o {params.block_tsv} -s {params.phased_stat}
        """


#### CONCAT SNPs ########
#########################

if config['gvcf_snv']:
    rule concatSNPs:
        """
        Rule to concat the identifed SNPs this will only be called by the user
        in case if he wanted to have only SNPs
        """
        input:
            vcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list),
            gvcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/data.{chr}.gvcf", aligner=wildcards.aligner, chr=chr_list),
        output:
            vcf = data_dir + "/snp/{aligner}/data.vcf",
            gvcf = data_dir + "/snp/{aligner}/data.gvcf",
        message: "Concat SNP files"
        benchmark: data_dir + "/benchmark/snp/{aligner}/concat_snp.txt"
        params:
            sample_name = SAMPLE_NAME,
        conda: VARIANT_ENV
        shell:"""
            echo "{params.sample_name}" > sample_name.txt && vcfcat {input.vcf} | vcfstreamsort | bcftools reheader --samples sample_name.txt -o {output.vcf} &&\
            vcfcat {input.gvcf} | vcfstreamsort | bcftools reheader --samples sample_name.txt -o {output.gvcf}
            """
else:
    rule concatSNPs:
        """
        Rule to concat the identified SNPs this will only be called by the user
        in case if he wanted to have only SNPs and indels and no gVCF required.
        """
        input:
            vcf = lambda wildcards: expand(data_dir + "/snp/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list),
        output:
            vcf = data_dir + "/snp/{aligner}/data.vcf"
        message: "Concat SNP files"
        benchmark: data_dir + "/benchmark/snp/{aligner}/concat_snp.txt"
        params:
            sample_name = SAMPLE_NAME,
        conda: VARIANT_ENV
        shell:"""
            echo "{params.sample_name}" > sample_name.txt && vcfcat {input.vcf} | vcfstreamsort | bcftools reheader --samples sample_name.txt -o {output.vcf}
            """

# #### Bgzip gVCF #########
# #########################

rule bgzipgVCFFile:
    """
    General rule to bgzip gVCF files
    """
    input:data_dir + "/{name}.gvcf"
    output:data_dir + "/{name}.gvcf.gz"
    threads: config['bgzip_threads']
    conda: VARIANT_ENV
    shell:"""
        bgzip -c -@ {threads} {input} > {output}
        """
