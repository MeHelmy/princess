
############################
######  PHASNIG RULES ######
###########################




#### GENOTYPING ####
####################

rule gt:
    """
    Genotype SNPs one chromosome per time.
    """
    input:
        bam=data_dir + "/align/{aligner}/data.bam",
        bam_index=data_dir + "/align/{aligner}/data.bam.bai",
        snps=data_dir + "/snp/{aligner}/data.{chr}.vcf",
    output:
        data_dir + "/gt/{aligner}/data.{chr}.vcf"
    params:
        reference=REFERENCES,
    conda: WHATSHAP_ENV
    log:
        data_dir + "/gt/{aligner}/data.{chr}.log"
    benchmark: data_dir + "/benchmark/gt/{aligner}/{chr}.benchmark.txt"
    shell:"""
        whatshap genotype --reference {params.reference} \
                      --ignore-read-groups \
                       --output {output} {input.snps} {input.bam} > {log} 2>&1
        """

#### PHASING ####
#################

rule phasing:
    """
    Phase SNPs one chromosome per time
    """
    input:
        bam=data_dir + "/align/{aligner}/data.bam",
        bam_index=data_dir + "/align/{aligner}/data.bam.bai",
        snps=data_dir + "/snp/{aligner}/data.{chr}.vcf",
    output:
        phased=temp(data_dir + "/phased/{aligner}/data.{chr}.vcf"),
    params:
        reference=REFERENCES,
        read_list=data_dir + "/phased/{aligner}/data.{chr}.reads",
    log:
        data_dir + "/phased/{aligner}/data.{chr}.log"
    conda: WHATSHAP_ENV
    benchmark: data_dir + "/benchmark/phase/{aligner}/{chr}.benchmark.txt"
    shell:"""
        whatshap phase --reference {params.reference} \
                     --output {output.phased} {input.snps} {input.bam} \
                     --ignore-read-groups \
                     --output-read-list {params.read_list} > {log} 2>&1
        """

#### CONCAT PHASING ####
########################

rule allPhased:
    """
    Concat all the phased SNPs into one file.
    """
    input:lambda wildcards: expand(data_dir + "/phased/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list),
    output: temp(data_dir + "/phased/{aligner}/data.vcf")
    conda: PRINCESS_ENV
    params:
        sample_name = SAMPLE_NAME,
    benchmark: data_dir + "/benchmark/phase/{aligner}/concat_phased.benchmark.txt"
    shell:"""
        echo "{params.sample_name}" > sample_name.txt && vcfcat {input} | vcfstreamsort | bcftools reheader --samples sample_name.txt -o {output}
        """

#### HAPLOTYPE BAM FILE ####
############################

rule partionBam:
    """
    Partion a bam file based on the phased SNPs,
    It will use the updated SNPs if the parental SNPs were provided.
    """
    input:
        bam = data_dir + "/align/{aligner}/data.bam",  # SM filed must be set to the sample name in vcf file
        bam_index = data_dir + "/align/{aligner}/data.bam.bai",
        snp = lambda wildcards: data_dir + "/phased/{aligner}/data_updated.vcf.gz" if config['update_snps'] else data_dir + "/phased/{aligner}/data.vcf.gz",
        snp_index = lambda wildcards: data_dir + "/phased/{aligner}/data_updated.vcf.gz.tbi" if config['update_snps'] else data_dir + "/phased/{aligner}/data.vcf.gz.tbi",
    output:
        hap_bam = data_dir + "/align/{aligner}/data_hap.bam"
    message: "Partioning bam file"
    conda: WHATSHAP_ENV
    params:
        ref = REFERENCES
    shell:"""
        whatshap haplotag -o {output.hap_bam} -r {params.ref} {input.snp} {input.bam}
        """
