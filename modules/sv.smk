
######################
###### SV RULES ######
#####################


#### SNIFFLES ####
##################
rule sniffles:
    """
    Identify structural variants using Sniffles2.
    """
    input:
        datain=data_dir + "/align/{aligner}/data_hap.bam" if config['phase_sv'] else data_dir + "/align/{aligner}/data.bam",
        data_index=data_dir + "/align/{aligner}/data_hap.bam.bai" if config['phase_sv'] else data_dir + "/align/{aligner}/data.bam.bai",
    output:
        dataout=data_dir + "/sv/{aligner}/sniffles.vcf",
        dataout_snf=data_dir + "/sv/{aligner}/sniffles.snf"
    message: "Running Sniffles in rule: {rule}\nUsing {input.datain} output:{output.dataout}"
    params:
        min_sv_len=config['min_sv_len'],
        sv_threads=config['sv_threads'],
        sample_name = SAMPLE_NAME,
        phase = "--phase" if config['phase_sv'] else "",
    conda: SNIFFLES_ENV
    priority: 2
    log: data_dir + "/sv/{aligner}/sniffles.log"
    benchmark: data_dir + "/benchmark/sv/{aligner}/sv.benchmark.txt"
    shell:"""
        sniffles --minsvlen {params.min_sv_len} --sample-id {params.sample_name} -t {params.sv_threads} --input {input.datain} --vcf {output.dataout} --snf {output.dataout_snf} {params.phase} > {log} 2>&1
        """

#### HAPLOTYPE SVs ####
#######################

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
        update_script = updat_sv,
    shell:"""
        python {params.update_script} {input.sv} {input.bam} {output} -c {params.min_conflict}
        """

#### SORTING SVs ####
#####################

rule vcfSort:
    """
    To concat the haplotype SVs with SNVs, SVs needs to be sorted first.
    """
    input:
        vcffile = data_dir + "/{sample}.vcf.gz",
        ref = REFERENCES,
    output:data_dir + "/{sample}.sorted.vcf.gz"
    conda: PRINCESS_ENV
    shell:"""
         zcat {input.vcffile} | awk 'BEGIN{{OFS="\t";}} /^#/{{print $0}} !/^#/{{ if ($2==0){{$2=1;print}} else {{print $0}} }}' |  bedtools sort -header -faidx {input.ref}.fai -i - | bgzip > {output}
        """
# bedtools sort -header -faidx {input.ref}.fai -i {input.vcffile} | bgzip > {output}
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
    conda: PRINCESS_ENV
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
        ref = REFERENCES,
    output: data_dir + "/sv/{aligner}/sv_snp.vcf.gz"
    params:
        extension = "z"  # output a compressed file.
    message: "Concat sv with SNPs"
    log: data_dir + "/sv/{aligner}/sv_snp.log"
    threads: config['samtools_threads']
    conda: PRINCESS_ENV
    shell:"""
        vcfcat  {input.sv} {input.snp}| bedtools sort -header -faidx {input.ref}.fai -i - | bgzip > {output} 2> {log}
        """
    # shell:"""
    #     bcftools concat -a -O {params.extension} -o {output} --threads {threads} {input.sv} {input.snp} > {log} 2>&1
    #     """
