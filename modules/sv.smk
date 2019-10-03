#####################
######  RULES ######
####################


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

#### SV STATISTICS ####
#######################

    input: expand("sv/{aligner}/sniffles.vcf", aligner=config['aligner'])
    output: "statitics/sv/data.stat"
    message: "calculating statistics for structural variant"
    benchmark: "benchmark/sv/stat.benchmark.txt"
    run:
        shell("""
        survivor stats {input} -1 -1 -1  {output}
        """)
