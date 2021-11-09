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
elif config['read_type'] == "ont":
    training_data=config["training_data_ont"]
elif config['read_type'] == "clr":
    training_data=config["training_data_clr"]
else:
    print("Unknow data type, supported format are: ont, ccs, and clr")
    exit(1)

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
rule callSNVsChunk:
    """
    Calling SNPs using clair
    """
    input:
        bam=data_dir + "/align/{aligner}/data.bam",
        data_index=data_dir + "/align/{aligner}/data.bam.bai",
        reference=REFERENCES,
    output:
        temp(data_dir + "/snp/{aligner}/{chrsplit}/chr.split.{chr}_{region,\d+}.vcf")
    params:
        train_data = training_data,
        minCoverage = config['clair_coverage'],
        start = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region)],
        end = lambda wildcards: chr_range[wildcards.chr][int(wildcards.region) + 1]
    benchmark: data_dir + "/benchmark/snp/{aligner}/{chrsplit}/{chr}_{region,\d+}.benchmark.txt"
    conda: CLAIR_ENV
    log: data_dir + "/snp/{aligner}/{chrsplit}/data.split.{chr}_{region}.log"
    # log: data_dir + "/snp/{aligner}/data.split.{chr,[A-Za-z0-9]+}_{region}.log"
    threads: config['clair_threads']
    shell:
        """
        export PATH=$PWD/bin/pypy/bin:$PATH && \
        clair.py callVarBam \
            --delay 0 \
            --chkpnt_fn {params.train_data} \
            --bam_fn {input.bam} \
            --ref_fn {input.reference} \
            --minCoverage {params.minCoverage} \
            --ctgName {wildcards.chr} \
            --ctgStart {params.start} \
            --ctgEnd {params.end} \
            --threads {threads} --call_fn {output}  > {log} 2>&1
        """


#### CALL VARINAT BY CHUNKS #######
###################################

rule concatChromosome:
    """
    Concat splited chromomsomes regions
    """
    input: lambda wildcards: expand(data_dir + "/snp/{aligner}/chrsplit/chr.split.{chr}_{region}.vcf", aligner=wildcards.aligner, chr=wildcards.chr, region=list(range(0,len(chr_range[wildcards.chr]) - 1))),
    output: temp(data_dir + "/snp/{aligner}/data.{chr}.vcf")
    message: "Concat variant split per Chromosome"
    conda: PRINCESS_ENV
    benchmark: data_dir + "/benchmark/snp/{aligner}/{chr}.benchmark.txt"
    params:
        temp_chr=data_dir + "/snp/{aligner}/data.{chr}_filtered.vcf",
        filter=config['filter_chrs'],
        read_type=config['read_type']
    shell:"""
        if [ {params.filter} == "True" ]; then
            find_max(){{
              filename=$(basename -- "$1")
              extension="${{1##*.}}"
              if [ "${{2}}" == "ont" ] || [ "${{2}}" == "clr" ]; then
                if [ "${{extension}}" == "gz" ]; then
                  zgrep -v "#" $1 | cut -f 6 | awk '$1 > 0 && $1 < 200 {{print}}' | sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -nr -k2,2 -k1,1 -t $'\t' | head -n1 |  awk '{{print $1}}'
                else
                  grep -v "#" $1 | cut -f 6 | awk '$1 > 0 && $1 < 200 {{print}}' | sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -nr -k2,2 -k1,1 -t $'\t' | head -n1 |  awk '{{print $1}}'
                fi
              elif [ "${{2}}" == "ccs" ]; then
                if [ "${{extension}}" == "gz" ]; then
                  zgrep -v "#" $1 | cut -f 6 | awk '$1 > 0 && $1 < 60 {{print}}' | sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -nr -k2,2 -k1,1 -t $'\t' | head -n1 |  awk '{{print $1}}'
                else
                  grep -v "#" $1 | cut -f 6 | awk '$1 > 0 && $1 < 60 {{print}}' | sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -nr -k2,2 -k1,1 -t $'\t' | head -n1 |  awk '{{print $1}}'
                fi
              else
                echo -e "Unknown technology ${{2}}"
              fi

            }}

            filsn () {{
            filename=$(basename -- "$1")
            extension="${{1##*.}}"
            if [[ -z ${{2:-}} ]];then
                min_qulaity=0
            else
              min_qulaity=$2
             fi
            if [ "${{extension}}" == "gz" ]; then
                zgrep -v "#" $1 |  cut -f 6  | awk -v min=$min_qulaity '$1 > min && $1 < 900 {{print}}'| sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -b -k2V -k1V | head -n1 | awk '{{print $1}}'
              else
                grep -v "#" $1 |  cut -f 6  | awk -v min=$min_qulaity '$1 > min && $1 < 900 {{print}}'| sort -n | uniq -c | awk '{{print $2,"\t",$1}}' | sort -b -k2V -k1V | head -n1 | awk '{{print $1}}'
            fi
            }}

            filecount=( {input} )
            count=${{#filecount[@]}}
            if [ "$count" -ge 2 ]; then
                vcfcat {input} | vcfstreamsort  > {params.temp_chr}\
                && first_max=$(find_max {params.temp_chr} {params.read_type})\
                && threshold=$(filsn {params.temp_chr} $first_max)\
                && awk -v threshold=$threshold '/^#/{{print}} !/^#/{{if ( $6 >= threshold ) {{print $0}}}}' {params.temp_chr} | awk '/^#/ {{ print }} !/^#/ {{ if ($4 != $5 ) {{ print }} }}' > {output}
            else
                if $(head -n 1000 {input} | grep -q -v "#") ; then
                    vcfcat {input}  | vcfstreamsort > {params.temp_chr}\
                    && first_max=$(find_max {params.temp_chr} {params.read_type})\
                    && threshold=$(filsn {params.temp_chr} $first_max)\
                    && awk -v threshold=$threshold '/^#/{{print}} !/^#/{{if ( $6 >= threshold ) {{print $0}}}}' {params.temp_chr} | awk '/^#/ {{ print }} !/^#/ {{ if ($4 != $5 ) {{ print }} }}' > {output}
                else
                    cp {input} {output}
                fi
            fi




        elif [ {params.filter} == "False" ]; then
            vcfcat {input} | awk '/^#/ {{ print }} !/^#/ {{ if ($4 != $5 ) {{ print }} }}' | vcfstreamsort | awk '/^#/ {{ print }} !/^#/ {{ if ($4 != $5 ) {{ print }} }}' > {output}
        else
            >&2 echo "Unknown option {params.filter}"
            exit 1
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
    message: "Indexing phacsed vcf file"
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
        phased_stat = data_dir + "/statistics/phased/phasing_stat.txt",
        block_tsv = data_dir + "/statistics/phased/blocks.tsv",
    benchmark: data_dir + "/benchmark/snp/{aligner}/update_snps.benchmark.txt"
    shell:"""
        mkdir -p statistics/phased  &&
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
