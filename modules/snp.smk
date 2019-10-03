#########################
######  SNPs RULES ######
#########################



# TODO: Clair needs a lot of configuration to run not just environement



#### CLAIR #######
##################

rule call_snps:
    """
    Calling SNPs using clair
    """
    input:
        bam="align/{aligner}/data.bam",
        data_index="align/{aligner}/data.bam.bai",
        reference=REFERENCES[ref[0]],
    output:
        "snp/{aligner}/data.{chr}.vcf"
    params:
        train_data=training_data,
        minCoverage=config['clair_coverage'],
        clair_location=CLAIR,
        process_script=process_clar,
        temp_out="{chr}_vcf.tmp",
        mypypy=config['clair_pypy'],
    benchmark: "benchmark/snp/{aligner}/{chr}.benchmark.txt"
    threads: config['clair_threads']
    shell:
        """
        export PATH=/users/mmahmoud/home/source/Clair/pypy2/pypy-7.0.0-linux_x86_64-portable/bin:$PATH && \
        source activate clar &&\
        python {params.clair_location} callVarBam \
            --chkpnt_fn {params.train_data} \
            --bam_fn {input.bam} \
            --ref_fn {input.reference} \
            --minCoverage {params.minCoverage} \
            --ctgName {wildcards.chr} \
            --threads {threads} --call_fn {output}
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
    input:"{sample}.vcf"
    output:"{sample}_update_header.vcf"
    message:"Update header file to change from float to string"
    run:
        shell("""
        sed 's/ID=PS,Number=1,Type=Integer,Descri/ID=PS,Number=1,Type=String,Descri/' {input} > {output}
        """)

#### BGZIP UPDATED FILE #######
###############################

# TODO:  The are a better way just by using {sample}.vcf as input , but debug needed cause it conflict with other rules.

rule bgzip_vcf:
    """
    bgzip phased file
    """
    input: "{sample}/data_update_header.vcf"
    output: "{sample}/data_update_header.vcf.gz"
    message: "bgzip the phased file"
    threads: config['bgzip_threads']
    run:
        shell("""
        bgzip -c -@ {threads} {input} > {output}
        """)

#### INDEXING VCF FILE ########
###############################

rule vcf_index:
    """
    Index VCF file.
    """
    input: "{sample}.vcf.gz"
    output: "{sample}.vcf.gz.tbi"
    message: "Indexing phacsed vcf file"
    conda:
    run:
        shell("""
        tabix -p vcf {input}
        """)

#### MERGING PHASED VCF FILE WITH PARENTAL SNPs ########
########################################################

rule merge_parental_snps:
    """
    If the user wanted to update identifed SNVs this will be the first rule in sequence,
    Input: phased SNVs after updating header using update_header the bgzip and index it using
    bgzip_vcf and vcf_index respectively.
    """
    input:
        sample_snps = "phased/{aligner}/data_update_header.vcf.gz",
        sample_snps_index = "phased/{aligner}/data_update_header.vcf.gz.tbi",
        maternal_snps = config['maternal_snps'],
        paternal_snps = config['paternal_snps'],
    output: "phased/{aligner}/data_paternal_maternal.vcf.gz"
    message: "merging vcf from samplepaternal and maternal respectively"
    benchmark: "benchmark/snp/{aligner}/merge_parental.benchmark.txt"
    conda:
    run:
        shell("""
        bcftools merge {input.sample_snps} {input.paternal_snps} {input.maternal_snps} | bgzip > {output}
        """)

#### UPDATING PHASED SNPs ########
##################################

rule update_snps:
    input: "phased/{aligner}/data_paternal_maternal.vcf"
    output:
        updated_vcf = "phased/{aligner}/data_updated.vcf",
    message: "Running update SNPs"
    params:
        update_script = config['updat_snps_script'],
        phased_stat = "statitics/phased/phasing_stat.txt",
        block_tsv = "statitics/phased/blocks.tsv",
    benchmark: "benchmark/snp/{aligner}/update_snps.benchmark.txt"
    conda:
    run:
        shell("""
        mkdir -p statitics/phased  &&
        python {params.update_script} -i {input} -u {output.updated_vcf} -o {params.block_tsv} -s {params.phased_stat}
        """)
