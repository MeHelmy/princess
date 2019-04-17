
rule gt:
    input:
        bam="align/{aligner}/data.bam",
        bam_index="align/{aligner}/data.bam.bai",
        snps="snp/{aligner}/data.{chr}.vcf",
    output:
        "gt/{aligner}/data.{chr}.vcf"
    params:
        reference=config['ref'],
    conda:
        WHATSHAP
    log:
        "gt/{aligner}/data.{chr}.log"
    benchmark:
        "gt/{aligner}/data.{chr}.benchmark.txt"
    shell:
        """
        whatshap genotype --reference {params.reference} \
                      --ignore-read-groups \
                       --output {output} {input.snps} {input.bam} > {log} 2>&1
        """

rule phasing:
    input:
        bam="align/{aligner}/data.bam",
        bam_index="align/{aligner}/data.bam.bai",
        snps="gt/{aligner}/data.{chr}.vcf",
    output:
        phased="phased/{aligner}/data.{chr}.vcf",
        read_list="phased/{aligner}/data.{chr}.reads",
    params:
        reference=config['ref'],
    conda:
        WHATSHAP
    log:
        "phased/{aligner}/data.{chr}.log"
    benchmark:
        "phased/{aligner}/data.{chr}.benchmark.txt"
    shell:
        """
        whatshap phase --reference {params.reference} \
                     --output {output.phased} {input.snps} {input.bam} \
                     --distrust-genotypes \
                     --ignore-read-groups \
                     --output-read-list {output.read_list} > {log} 2>&1
        """

rule all_phased:
    input:lambda wildcards: expand("phased/{aligner}/data.{chr}.vcf", aligner=wildcards.aligner, chr=chr_list[ref[0]]),
    output:
        "phased/{aligner}/data.vcf"
    conda:
        WHATSHAP
    shell:
        "vcfcat {input} | vcfstreamsort > {output}"
