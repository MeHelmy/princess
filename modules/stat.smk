#####################
######  RULES ######
####################


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
