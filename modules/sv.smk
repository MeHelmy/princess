######################
###### SV RULES ######
#####################


#### SNIFFLES ####
##################

rule sniffles:
    """
    Identify structural variants using Sniffles.
    """
    input:
        datain="align/{aligner}/data.bam",
        data_index="align/{aligner}/data.bam.bai",
    output:
        dataout="sv/{aligner}/sniffles.vcf"
    message: "Running Sniffles"
    params:
        coverage=config['sniffles_coverage'],
    benchmark: "sv/{aligner}/sniffles.benchmark.txt"
    conda:
        SV
    log: "sv/{aligner}/sniffles.log"
    benchmark: "benchmark/sv/{aligner}/sv.benchmark.txt"
    run:
        shell("""
        sniffles --min_support {params.coverage} --mapped_reads {input.datain} --vcf {output.dataout} --num_reads_report -1 --genotype --report_seq  > {log} 2>&1
        """)

#### HAPLOTYPE SVs ####
#######################

# TODO: rule is not finished script needs to be updated to change hp filed according to whashap results
rule phase_sv:
    """
    This rules takes as input a taged tabed bam file from whatshap and vcf file contains SVs and update the SVs to add haplotype HP and phase blocks PS.
    """
    input:
        bam = "align/{aligner}/data_hap.tab",
        sv = "sv/{aligner}/sniffles.vcf",
    output:"sv/{aligner}/sniffles_hp_updated.vcf"
    message: "Updating SVs using align/{aligner}/data_hap.tab"
    params:
        update_script = updat_sv
    run:
        shell("""
        python {params.update_script} {input.sv} {input.bam} {output}
        """)

#### SORTING SVs ####
#####################

rule vcf_sort:
    """
    To concat the haplotype SVs with SNVs, SVs needs to be sorted first.
    """
    input:"{sample}.vcf.gz"
    output:"{sample}.sorted.vcf.gz"
    conda: ALIGN
    run:
        shell("""
        bcftools sort -O z -o {output} {input}
        """)

#### BGZIP SVs ####
###################

rule bgzip_file:
    """
    General rule to bgzip files
    """
    input:"{name}.vcf"
    output:"{name}.vcf.gz"
    threads: config['bgzip_threads']
    run:
        shell("""
        bgzip -c -@ {threads} {input} > {output}
        """)

#### CHANGE SVs SAMPLE NAME ####
###############################

rule change_sample_name:
    """
    Sniffles name the sample as the bam, but Clair call it SAMPLE
    This rule will change the sample name in the SV file.
    """
    input:"{sample}.sorted.vcf.gz"
    output:"{sample}.sorted.namechnage.vcf.gz"
    conda: ALIGN
    run:
        shell("""
        echo SAMPLE > sample.name && bcftools reheader -s sample.name  -o {output} {input} && rm sample.name
        """)

#### CONCAT SVs WITH SNPs ####
##############################

rule SV_SNP_compained:
    """
    Concat haplotyped SNPs with haplotyped and Genotyped SVs.
    """
    input:
        sv = "sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz",
        sv_index ="sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz.tbi",
        snp =  lambda wildcards: "phased/{aligner}/data_updated.vcf.gz" if config['update_snps'] else "phased/{aligner}/data.vcf.gz",
    output: "sv/{aligner}/sv_snp.vcf.gz"
    params:
        extension = "z"  # output a compressed file.
    message: "Concat sv with SNPs"
    log: "sv/{aligner}/sv_snp.log"
    threads: config['samtools_threads']
    conda: ALIGN
    run:
        shell("""
        bcftools concat -a -O {params.extension} -o {output} --threads {threads} {input.sv} {input.snp} > {log} 2>&1
        """)
