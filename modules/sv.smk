
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
        datain=data_dir + "/align/{aligner}/data.bam",
        data_index=data_dir + "/align/{aligner}/data.bam.bai",
    output:
        dataout=data_dir + "/sv/{aligner}/sniffles.vcf"
    message: "Running Sniffles"
    params:
        coverage=config['sniffles_coverage'],
    # benchmark: data_dir + "/sv/{aligner}/sniffles.benchmark.txt"
    conda: PRINCESS_ENV
    priority: 1
    log: data_dir + "/sv/{aligner}/sniffles.log"
    benchmark: data_dir + "/benchmark/sv/{aligner}/sv.benchmark.txt"
    shell:"""
        sniffles --min_support {params.coverage} --mapped_reads {input.datain} --vcf {output.dataout} --num_reads_report -1 --genotype  > {log} 2>&1
        """

#### HAPLOTYPE SVs ####
#######################

# TODO: rule is not finished script needs to be updated to change hp filed according to whashap results
rule phaseSVs:
    """
    This rules takes as input a taged tabed bam file
    from whatshap and vcf file contains SVs and update
    the SVs to add haplotype HP and phase blocks PS.
    """
    input:
        bam = data_dir + "/align/{aligner}/data_hap.tab",
        sv = data_dir + "/sv/{aligner}/sniffles.vcf",
    output:data_dir + "/sv/{aligner}/sniffles_hp_updated.vcf"
    message: "Updating SVs using align/{aligner}/data_hap.tab"
    params:
        update_script = updat_sv
    shell:"""
        python {params.update_script} {input.sv} {input.bam} {output}
        """

#### SORTING SVs ####
#####################

rule vcfSort:
    """
    To concat the haplotype SVs with SNVs, SVs needs to be sorted first.
    """
    input:data_dir + "/{sample}.vcf.gz"
    output:data_dir + "/{sample}.sorted.vcf.gz"
    conda: PRINCESS_ENV
    shell:"""
        zcat {input} | vcfsort | bgzip > {output}
        """
    # shell:"""
    #     bcftools sort -O z -o {output} {input}
    #     """

#### BGZIP SVs ####
###################

rule bgzipFile:
    """
    General rule to bgzip files
    """
    input:data_dir + "/{name}.vcf"
    output:data_dir + "/{name}.vcf.gz"
    threads: config['bgzip_threads']
    shell:"""
        bgzip -c -@ {threads} {input} > {output}
        """

#### CHANGE SVs SAMPLE NAME ####
###############################

rule changeSampleName:
    """
    Sniffles name the sample as the bam, but Clair call it SAMPLE
    This rule will change the sample name in the SV file.
    """
    input:data_dir + "/{sample}.sorted.vcf.gz"
    output:data_dir + "/{sample}.sorted.namechnage.vcf.gz"
    conda: PRINCESS_ENV
    shell:"""
        echo SAMPLE > sample.name && bcftools reheader -s sample.name  -o {output} {input} && rm sample.name
        """

#### CONCAT SVs WITH SNPs ####
##############################

rule SVsSNPsCompained:
    """
    Concat haplotyped SNPs with haplotyped and Genotyped SVs.
    """
    input:
        sv = data_dir + "/sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz",
        sv_index =data_dir + "/sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz.tbi",
        snp =  lambda wildcards: data_dir + "/phased/{aligner}/data_updated.vcf.gz" if config['update_snps'] else data_dir + "/phased/{aligner}/data.vcf.gz",
    output: data_dir + "/sv/{aligner}/sv_snp.vcf.gz"
    params:
        extension = "z"  # output a compressed file.
    message: "Concat sv with SNPs"
    log: data_dir + "/sv/{aligner}/sv_snp.log"
    threads: config['samtools_threads']
    conda: PRINCESS_ENV
    shell:"""
        bcftools concat -a -O {params.extension} -o {output} --threads {threads} {input.sv} {input.snp} > {log} 2>&1
        """
