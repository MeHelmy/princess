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
        bam="align/{aligner}/data.bam",
        bam_index="align/{aligner}/data.bam.bai",
        snps="snp/{aligner}/data.{chr}.vcf",
    output:
        temp("gt/{aligner}/data.{chr}.vcf")
    params:
        reference=REFERENCES[ref[0]],
    conda:
        WHATSHAP
    log:
        "gt/{aligner}/data.{chr}.log"
    benchmark: "benchmark/gt/{aligner}/{chr}.benchmark.txt"
    run:
        shell("""
        whatshap genotype --reference {params.reference} \
                      --ignore-read-groups \
                       --output {output} {input.snps} {input.bam} > {log} 2>&1
        """)

#### PHASING ####
#################

rule phasing:
    """
    Phase SNPs one chromosome per time
    """
    input:
        bam="align/{aligner}/data.bam",
        bam_index="align/{aligner}/data.bam.bai",
        snps="gt/{aligner}/data.{chr}.vcf",
    output:
        phased=temp("phased/{aligner}/data.{chr}.vcf"),
    params:
        reference=REFERENCES[ref[0]],
        read_list="phased/{aligner}/data.{chr}.reads",
    log:
        "phased/{aligner}/data.{chr}.log"
    conda:
        WHATSHAP
    benchmark: "benchmark/phase/{aligner}/{chr}.benchmark.txt"
    run:
        shell("""
        whatshap phase --reference {params.reference} \
                     --output {output.phased} {input.snps} {input.bam} \
                     --distrust-genotypes \
                     --ignore-read-groups \
                     --output-read-list {params.read_list} > {log} 2>&1
        """)

#### CONCAT PHASING ####
########################

rule all_phased:
    """
    Concat all the phased SNPs into one file.
    """
    input:lambda wildcards: expand("phased/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list[ref[0]]),
    output: "phased/{aligner}/data.vcf"
    conda:
        WHATSHAP
    benchmark: "benchmark/phase/{aligner}/concat_phased.benchmark.txt"
    run:
        shell("""
        vcfcat {input} | vcfstreamsort > {output}
        """)

#### HAPLOTYPE BAM FILE ####
############################

rule partion_bam:
    """
    Partion a bam file based on the phased SNPs,
    It will use the updated SNPs if the parental SNPs were provided.
    """
    input:
        bam = "align/{aligner}/data.bam",  # SM filed myst be set to the sample name in vcf file
        bam_index = "align/{aligner}/data.bam.bai",
        snp = lambda wildcards: "phased/{aligner}/data_updated.vcf.gz" if config['update_snps'] else "phased/{aligner}/data.vcf.gz",
        snp_index = lambda wildcards: "phased/{aligner}/data_updated.vcf.gz.tbi" if config['update_snps'] else "phased/{aligner}/data.vcf.gz.tbi",
        #"phased/{aligner}/data_updated.vcf.gz" if config['update_snps'] else "phased/{aligner}/data.vcf.gz
    output:
        hap_bam = "align/{aligner}/data_hap.bam"
    message: "Partioning bam file"
    params:
        ref = REFERENCES[ref[0]]
    run:
        shell("""
        whatshap haplotag -o {output.hap_bam} -r {params.ref} {input.snp} {input.bam}
        """)
