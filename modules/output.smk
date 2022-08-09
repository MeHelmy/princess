
##########################
######  OUTPUT RULES ####
#########################

#### Align Moving ########
#########################

rule mvAlign:
    input:
        bam = data_dir + "/align/{aligner}/data.bam",
        bamindex = data_dir + "/align/{aligner}/data.bam.bai",
    output:
        bam = data_dir +'/result' + '/.aligning.{aligner}.done',
    message: "Moving Aligned bam to result directory {input.bam}"
    params:
        bam = data_dir +'/result' + "/aligning.{sample}.{{aligner}}.bam".format(sample=SAMPLE_NAME),
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2
    }}
    mv {input.bamindex} {params.bam}.bai && mv {input.bam} {params.bam} &&\
    mkdir -p {data_dir}/log &&\
    mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || : &&\
    rm_last2 {input.bam} &&\
    touch {output}
    """


#### SVs Moving ########
########################

rule mvSV:
    input:
        vcf = data_dir + "/sv/{aligner}/sniffles.vcf",
        snf = data_dir + "/sv/{aligner}/sniffles.snf",
        bam = data_dir + "/align/{aligner}/data.bam",
        bamindex = data_dir + "/align/{aligner}/data.bam.bai",
    output:
        vcf = data_dir +'/result' + '/.SVs.{aligner}.done',
    params:
        vcf = data_dir +'/result' + '/SVs.{aligner}.vcf',
        snf = data_dir +'/result' + '/SVs.{aligner}.snf',
        bam = data_dir +'/result' + '/aligning.{aligner}.bam',
        bamindex = data_dir +'/result' + '/aligning.{aligner}.bam.bai',
    message: "Moving called SVs to result directory {input}"
    priority: 1
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2
    }}
    mkdir -p {data_dir}/log &&\
    if [ -f {data_dir}/sv/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/sv/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    mv {input.vcf} {params.vcf} && rm_last2 {input.vcf} &&\
    mv {input.snf} {params.snf} && rm_last2 {input.snf} &&\
    mv {input.bam} {params.bam} && mv {input.bamindex} {params.bamindex} && rm_last2 {input.bam} &&\
    touch {output.vcf}
    """


#### SNVs Moving ########
#########################

rule mvSNV:
    input:
        vcf = data_dir + "/snp/{aligner}/data.vcf.gz",
        vcfindex = data_dir + "/snp/{aligner}/data.vcf.gz.tbi",
        bam = data_dir + "/align/{aligner}/data.bam",
        bamindex = data_dir + "/align/{aligner}/data.bam.bai",
    output:
        data_dir +'/result' + '/.SNVs.{aligner}.done'
    params:
        vcf = data_dir +'/result' + '/SNVs.{aligner}.vcf.gz',
        vcfindex = data_dir +'/result' + '/SNVs.{aligner}.vcf.gz.tbi',
        bam = data_dir +'/result' + '/aligning.{aligner}.bam',
        bamindex = data_dir +'/result' + '/aligning.{aligner}.bam.bai',
    message: "Moving called SNVs to result directory {input}"
    priority: 1
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2
    }}
    mkdir -p {data_dir}/log &&\
    mv {input.vcf} {params.vcf} &&\
    mv {input.vcfindex} {params.vcfindex}  &&\
    mv {input.bam} {params.bam} &&\
    mv {input.bamindex} {params.bamindex}  &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    rm_last2 {input.vcf} && rm_last2 {input.bam} &&\
    touch {output}
    """

#### Variants Moving ########
############################

rule mvVariants:
    input:
        snv = data_dir + "/snp/{aligner}/data.vcf.gz",
        snvindex = data_dir + "/snp/{aligner}/data.vcf.gz.tbi",
        sv = data_dir + "/sv/{aligner}/sniffles.vcf",
        snf = data_dir + "/sv/{aligner}/sniffles.snf",
        bam = data_dir + "/align/{aligner}/data.bam",
        bamindex = data_dir + "/align/{aligner}/data.bam.bai",
    output:
        data_dir + "/result" + "/.varaint.{aligner}.done"
    message: "Moving called SNVs to result directory {input}"
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2
    }}
    mkdir -p {data_dir}/log &&\
    if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/sv/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/sv/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/chrsplit/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/chrsplit/*.log {data_dir}/log || :
    fi &&\
    mv {input.snv} {data_dir}/result/SNVs.{wildcards.aligner}.vcf.gz &&\
    mv {input.snvindex} {data_dir}/result/SNVs.{wildcards.aligner}.vcf.gz.tbi &&\
    rm_last2 {input.snv} &&\
    mv {input.sv} {data_dir}/result/SVs.{wildcards.aligner}.vcf &&\
    mv {input.snf} {data_dir}/result/SVs.{wildcards.aligner}.snf &&\
    rm_last2 {input.sv} &&\
    rm_last2 {input.snf} &&\
    mv {input.bam} {data_dir}/result/align.{wildcards.aligner}.bam &&\
    mv {input.bamindex} {data_dir}/result/align.{wildcards.aligner}.bam.bai &&\
    rm_last2 {input.bam} &&\
    touch {output}
    """

#### Phasing Moving ########
############################

rule mvPhasing:
    input:
        snv = data_dir + "/phased/{aligner}/data.vcf.gz",
        snvindex = data_dir + "/phased/{aligner}/data.vcf.gz.tbi",
        bam = data_dir + "/align/{aligner}/data.bam",
        bamindex = data_dir + "/align/{aligner}/data.bam.bai",
    output:
        snv = data_dir + "/result" + "/phased.SNVs.{aligner}.done",
    params:
        snv = data_dir + "/result" + "/phased.SNVs.{aligner}.vcf.gz",
        snvindex = data_dir + "/result" + "/phased.SNVs.{aligner}.vcf.gz.tbi",
        bam = data_dir +'/result' + '/aligning.{aligner}.bam',
        bamindex = data_dir +'/result' + '/aligning.{aligner}.bam.bai',
    message: "Moving called phased SNVs to result directory {input}"
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2 ||:
    }}
    function rm_last1() {{
    d1=$(dirname $1)
    rm -rf $d1 ||:
    }}
    mkdir -p {data_dir}/log &&\
    if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/phased/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/phased/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/phased/{wildcards.aligner}/*.txt ]; then
        mv {data_dir}/phased/{wildcards.aligner}/*.txt {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/chrsplit/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/chrsplit/*.log {data_dir}/log || :
    fi &&\
    mv {input.snv} {params.snv} &&\
    mv {input.snvindex} {params.snvindex} &&\
    mv {input.bam} {params.bam} &&
    mv {input.bamindex} {params.bamindex} &&\
    rm_last2 {input.snv} &&\
    rm_last2 {input.bam} &&\
    rm_last1 {data_dir}/snp/{wildcards.aligner} &&\
    touch {output}
    """

#### All Moving ########
########################

rule mvmethylation:
    input:
        methylation = data_dir + "/meth/"+ "{aligner}" + "/methylation_calls_hap.tsv",
    output:
        methylation = data_dir + "/result" + "/methylation.{aligner}_calls_hap.tsv"
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2 ||:
    }}
    mv {input.methylation} {output.methylation} &&\
    rm_last2 input.methylation}
    """

rule mvParentalPhased:
    input:
        stat = data_dir + "/stat.txt" if config['sample_list'] else data_dir + "/stat.NoReads.txt",
        phasedSNVs = data_dir + "/phased/{aligner}/data_updated.vcf",
        # phasedSVs = data_dir + "/sv/{aligner}/sniffles_hp_updated.vcf",
        bam = data_dir + "/align/{aligner}/data_hap.bam",
        bamindex = data_dir + "/align/{aligner}/data_hap.bam.bai",
    output:
        stat = data_dir + "/result/.allReadsparental.{aligner}.txt"  #if config['sample_list'] else data_dir + "/result/.allNoReadsparental.{aligner}.txt",
    params:
        stat = data_dir + "/result/stat.{aligner}.txt",
        phasedSNVs = data_dir + "/result/{aligner}.phased.SNVs.vcf",
        # phasedSVs = data_dir + "/result/{aligner}.phased.SVs.vcf",
        bam = data_dir + "/result/{aligner}.hap.bam",
        bamindex = data_dir + "/result/{aligner}.hap.bam.bai",
    shell:"""
    function rm_last2() {{
    d1=$(dirname $1)
    d2=$(dirname $d1)
    rm -rf $d2 ||:
    }}
    function rm_last1() {{
    d1=$(dirname $1)
    rm -rf $d1 ||:
    }}
    mkdir -p {data_dir}/log &&\
    if [ -f {data_dir}/sv/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/sv/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    if [ -f {data_dir}/phased/{wildcards.aligner}/*.log ]; then
        mv {data_dir}/phased/{wildcards.aligner}/*.log {data_dir}/log || :
    fi &&\
    mv {input.stat} {params.stat} &&\
    mv {data_dir}/statitics {data_dir}/result &&\
    mv {input.phasedSNVs} {params.phasedSNVs} &&\
    bgzip {params.phasedSNVs} &&\
    tabix {params.phasedSNVs}.gz &&\
    mv {input.bam} {params.bam} &&\
    mv {input.bamindex} {params.bamindex} &&\
    rm_last2 {input.phasedSNVs} &&\
    rm_last2 {input.bam} &&\
    touch {output}
    """
# mv {input.phasedSVs} {params.phasedSVs} &&\
# rm_last2 {input.phasedSVs} &&\

if config['gvcf_snv']:
    rule mvNoParentalPhased:
        """
        The rule here will make sure to get the rule for all command without parental comparison and later will delete the data.
        """
        input:
            stat = data_dir + "/stat.txt" if config['sample_list'] else data_dir + "/stat.NoReads.txt",
            # phasedSvsSNVs = data_dir + "/sv/{aligner}/sv_snp.vcf.gz",
            phasedSNVs = data_dir + "/phased/{aligner}/data.vcf.gz",
            phasedSNVsindex = data_dir + "/phased/{aligner}/data.vcf.gz.tbi",
            SVs = data_dir + "/sv/{aligner}/sniffles.vcf",
            snf = data_dir + "/sv/{aligner}/sniffles.snf",
            # phasedSVs = data_dir + "/sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz",
            bam = data_dir + "/align/{aligner}/data_hap.bam",
            bamindex = data_dir + "/align/{aligner}/data_hap.bam.bai",
            gvcf = data_dir + "/snp/{aligner}/data.gvcf.gz",
        output:
            stat = data_dir + "/result/.all.Reads.{aligner}.txt" if config['sample_list'] else data_dir + "/result/.all.noReads.{aligner}.txt",
        params:
            stat = data_dir + "/result/stat.{sample}.{{aligner}}.txt".format(sample=SAMPLE_NAME),
            # phasedSvsSNVs = data_dir + "/result/{aligner}.phased.sv_snp.vcf.gz",
            phasedSNVs = data_dir + "/result/{sample}.{{aligner}}.phased.SNVs.vcf.gz".format(sample=SAMPLE_NAME),
            phasedSNVsindex = data_dir + "/result/{sample}.{{aligner}}.phased.SNVs.vcf.gz.tbi".format(sample=SAMPLE_NAME),
            SVs = data_dir + "/result/{sample}.{{aligner}}.SVs.vcf".format(sample=SAMPLE_NAME),
            snf = data_dir + "/result/{sample}.{{aligner}}.SVs.snf".format(sample=SAMPLE_NAME),
            # phasedSVs = data_dir + "/result/{aligner}.SVs.phased.vcf.gz",
            bam = data_dir + "/result/{sample}.{{aligner}}.hap.bam".format(sample=SAMPLE_NAME),
            bamindex = data_dir + "/result/{sample}.{{aligner}}.hap.bam.bai".format(sample=SAMPLE_NAME),
            copy_gvcf = "True" if config['gvcf_snv'] else "False",
            gvcf = data_dir + "/result/{sample}.{{aligner}}.SNVs.gvcf.gz".format(sample=SAMPLE_NAME),
        shell:"""
        function rm_last2() {{
        d1=$(dirname $1)
        d2=$(dirname $d1)
        rm -rf $d2 ||:
        }}
        function rm_last1() {{
        d1=$(dirname $1)
        rm -rf $d1 ||:
        }}
        mkdir -p {data_dir}/log &&\
        if [ -f {data_dir}/sv/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/sv/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/phased/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/phased/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        mv {input.gvcf} {params.gvcf} &&\
        mv {input.stat} {params.stat} &&\
        mv {input.phasedSNVs} {params.phasedSNVs} &&\
        mv {input.phasedSNVsindex} {params.phasedSNVsindex} &&\
        mv {input.SVs} {params.SVs} &&\
        mv {input.snf} {params.snf} &&\
        mv {input.bam} {params.bam} &&\
        mv {input.bamindex} {params.bamindex} &&\
        rm_last2 {input.SVs} &&\
        rm_last2 {input.phasedSNVs} &&\
        rm_last2 {input.bam} &&\
        rm_last1 {data_dir}/snp/{wildcards.aligner} &&\
        touch {output}
        """
else:
    rule mvNoParentalPhased:
        """
        The rule here will make sure to get the rule for all command without parental comparison and later will delete the data.
        """
        input:
            stat = data_dir + "/stat.txt" if config['sample_list'] else data_dir + "/stat.NoReads.txt",
            # phasedSvsSNVs = data_dir + "/sv/{aligner}/sv_snp.vcf.gz",
            phasedSNVs = data_dir + "/phased/{aligner}/data.vcf.gz",
            phasedSNVsindex = data_dir + "/phased/{aligner}/data.vcf.gz.tbi",
            SVs = data_dir + "/sv/{aligner}/sniffles.vcf",
            snf = data_dir + "/sv/{aligner}/sniffles.snf",
            # phasedSVs = data_dir + "/sv/{aligner}/sniffles_hp_updated.sorted.namechnage.vcf.gz",
            bam = data_dir + "/align/{aligner}/data_hap.bam",
            bamindex = data_dir + "/align/{aligner}/data_hap.bam.bai",
        output:
            stat = data_dir + "/result/.all.Reads.{aligner}.txt" if config['sample_list'] else data_dir + "/result/.all.noReads.{aligner}.txt",
        params:
            stat = data_dir + "/result/stat.{sample}.{{aligner}}.txt".format(sample=SAMPLE_NAME),
            # phasedSvsSNVs = data_dir + "/result/{aligner}.phased.sv_snp.vcf.gz",
            phasedSNVs = data_dir + "/result/{sample}.{{aligner}}.phased.SNVs.vcf.gz".format(sample=SAMPLE_NAME),
            phasedSNVsindex = data_dir + "/result/{sample}.{{aligner}}.phased.SNVs.vcf.gz.tbi".format(sample=SAMPLE_NAME),
            SVs = data_dir + "/result/{sample}.{{aligner}}.SVs.vcf".format(sample=SAMPLE_NAME),
            snf = data_dir + "/result/{sample}.{{aligner}}.SVs.snf".format(sample=SAMPLE_NAME),
            # phasedSVs = data_dir + "/result/{aligner}.SVs.phased.vcf.gz",
            bam = data_dir + "/result/{sample}.{{aligner}}.hap.bam".format(sample=SAMPLE_NAME),
            bamindex = data_dir + "/result/{sample}.{{aligner}}.hap.bam.bai".format(sample=SAMPLE_NAME),
            copy_gvcf = "True" if config['gvcf_snv'] else "False",
        shell:"""
        function rm_last2() {{
        d1=$(dirname $1)
        d2=$(dirname $d1)
        rm -rf $d2 ||:
        }}
        function rm_last1() {{
        d1=$(dirname $1)
        rm -rf $d1 ||:
        }}
        mkdir -p {data_dir}/log &&\
        if [ -f {data_dir}/sv/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/sv/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/snp/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/snp/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/align/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/align/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        if [ -f {data_dir}/phased/{wildcards.aligner}/*.log ]; then
            mv {data_dir}/phased/{wildcards.aligner}/*.log {data_dir}/log || :
        fi &&\
        mv {input.stat} {params.stat} &&\
        mv {input.phasedSNVs} {params.phasedSNVs} &&\
        mv {input.phasedSNVsindex} {params.phasedSNVsindex} &&\
        mv {input.SVs} {params.SVs} &&\
        mv {input.snf} {params.snf} &&\
        mv {input.bam} {params.bam} &&\
        mv {input.bamindex} {params.bamindex} &&\
        rm_last2 {input.SVs} &&\
        rm_last2 {input.phasedSNVs} &&\
        rm_last2 {input.bam} &&\
        rm_last1 {data_dir}/snp/{wildcards.aligner} &&\
        touch {output}
        """
# mv {input.phasedSVs} {params.phasedSVs} &&\
# mv {input.phasedSvsSNVs} {params.phasedSvsSNVs} &&\
# rm_last2 {input.phasedSvsSNVs} &&\
