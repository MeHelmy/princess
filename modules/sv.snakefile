rule sniffles:
    """
    Identify structure variation using Sniffles.
    """
    input:
        datain="align/{aligner}/data.bam",
        data_index="align/{aligner}/data.bam.bai",
    output:
        dataout="sv/{aligner}/sniffles.vcf"
    message:
        "Running Sniffles"
    params:
        coverage=config['sniffles_coverage'],
    benchmark:
        "sv/{aligner}/sniffles.benchmark.txt"
    conda:
        SV
    log:
        "sv/{aligner}/sniffles.log"
    run:
        shell("""
            sniffles --min_support {params.coverage} --mapped_reads {input.datain} --vcf {output.dataout} --num_reads_report -1 --genotype --report_seq  > {log} 2>&1
            """)
