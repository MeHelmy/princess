###############################
######  STATISTICS RULES ######
###############################



#### RAW READS STATISTICS ####
##############################

rule reads_stat:
    """
    Input is the reads in directory output is info about reads
    """
    input: expand(data_dir + "/{sample}", sample=sample_list)
    output:data_dir + "/statitics/raw_reads/reads_stat.txt",
    message: "Calculating read coverage statitics for: {input}",
    params:
        read_stat_script = rawcoverage_script,
    threads: config['read_raw_coverage_threads']
    benchmark: data_dir + "/benchmark/raw_reads/stat.benchmark.txt"
    conda: PRINCESS_ENV
    shell:
        """
        python {params.read_stat_script} -i {input} -o {output} -t {threads}
        """

#### BAM STATISTICS ####
########################

rule bam_stat:
    """
    Calculate statistics from merged bam file
    """
    input:data_dir + "/align/{aligner}/data.bam"
    output:data_dir + "/statitics/{aligner}/data.stat"
    message:"Calculating aligned reads statistics from bam file"
    benchmark: data_dir + "/benchmark/align/{aligner}/stat.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        samtools stats {input} > {output}
        """

#### SV STATISTICS ####
#######################

rule sv_stat:
    input: expand(data_dir + "/sv/{aligner}/sniffles.vcf", aligner=config['aligner'])
    output: data_dir + "/statitics/sv/data.stat"
    message: "calculating statistics for structural variant"
    benchmark: data_dir + "/benchmark/sv/stat.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        SURVIVOR stats {input} -1 -1 -1  {output}
        """

#### SNPs STATISTICS ####
#########################

rule snp_stat:
    input:
        snp_file = expand(data_dir + "/phased/{aligner}/data.vcf.gz", aligner=config['aligner']) ,
        snp_file_index = expand(data_dir + "/phased/{aligner}/data.vcf.gz.tbi", aligner=config['aligner']) ,
    output: data_dir + "/statitics/snp/snp.txt",
    message: "Calculate SNPs statistics"
    benchmark: data_dir + "/benchmark/snp/stat.benchmark.txt"
    conda: PRINCESS_ENV
    shell:"""
        bcftools stats {input.snp_file} > {output}
        """


#### ALL STATISTICS ####
#######################

rule stat:
    input:
        expand(data_dir + "/statitics/{aligner}/data.stat", aligner=config['aligner']),
        data_dir + "/statitics/raw_reads/reads_stat.txt",\
        data_dir + "/statitics/sv/data.stat",\
        data_dir + "/statitics/snp/snp.txt",
    output: data_dir + "/stat.txt"
    shell: "touch {output}"
