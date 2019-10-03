###############################
######  STATISTICS RULES ######
###############################



#### RAW READS STATISTICS ####
##############################

rule reads_stat:
    """
    Input is the reads in directory output is info about reads
    """
    input: expand("{sample}", sample=sample_list)
    output:"statitics/raw_reads/reads_stat.txt",
    message: "Calculating read coverage statitics for: {input}",
    params:
        read_stat_script = rawcoverage_script,
    threads: config['read_raw_coverage_threads']
    benchmark: "benchmark/raw_reads/stat.benchmark.txt"
    run:
        shell(
        """
        python {params.read_stat_script} -i {input} -o {output} -t {threads}
        """
        )

#### BAM STATISTICS ####
########################

rule bam_stat:
    """
    Calculate statistics from merged bam file
    """
    input:"align/{aligner}/data.bam"
    output:"statitics/{aligner}/data.stat"
    message:"Calculating aligned reads statistics from bam file"
    benchmark: "benchmark/align/{aligner}/stat.benchmark.txt"
    conda: ALIGN
    run:
        shell("""
        samtools stats {input} > {output}
        """)

#### SV STATISTICS ####
#######################

rule sv_stat:
    input: expand("sv/{aligner}/sniffles.vcf", aligner=config['aligner'])
    output: "statitics/sv/data.stat"
    message: "calculating statistics for structural variant"
    benchmark: "benchmark/sv/stat.benchmark.txt"
    run:
        shell("""
        survivor stats {input} -1 -1 -1  {output}
        """)

#### SNPs STATISTICS ####
#########################

rule snp_stat:
    input:
        snp_file = expand("phased/{aligner}/data.vcf.gz", aligner=config['aligner']) ,
        snp_file_index = expand("phased/{aligner}/data.vcf.gz.tbi", aligner=config['aligner']) ,
    output: "statitics/snp/snp.txt",
    message: "Calculate SNPs statistics"
    benchmark: "benchmark/snp/stat.benchmark.txt"
    run:
        shell("""
        bcftools stats {input.snp_file} > {output}
        """)


#### ALL STATISTICS ####
#######################

rule stat:
    input:
        expand("statitics/{aligner}/data.stat", aligner=config['aligner']),
        "statitics/raw_reads/reads_stat.txt",\
        "statitics/sv/data.stat",\
        "statitics/snp/snp.txt",
    output: "stat.txt"
    shell: "touch {output}"
