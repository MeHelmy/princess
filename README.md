## What is new?
- Clair3 for calling single nucleotide polymorphisms (SNPs) and insertions/deletions (indels)
  - Ability to use different models than the default one that comes with Clair3, which can be helpful in cases where there is new kit/training dataset or when working with data other than the human genome.
- Sniffles2 for detecting structural variants (SVs)
- Generation of a gVCF file for cohort analysis
- Generation of an SNF file for cohort structural variant analysis
- The pipeline has been fully tested on both PBS and Slurm systems with easy configuration
- The main conda environment has been updated for improved granularity.
---

Princess is a fast and scalable framework to detect and report haplotype resolved Single Nucleotide Variants (SNV) and Structural Variations (SVs) at scale. It can leverage your cluster environment to speed up the detection which starts with one or many fasta or fastq files.

![princess](./pictures/leia.jpg)

## Princess

* __Mapping__:  NGMLR or Minimap2
* __SNVs__: Clair3 
* __SVs__: Sniffles2
* __Phasing SNVs__: WhatsHap
* __Phasing SVs__: PRINCESS-subtool
* __Extend Phasing__: PRINCESS-subtool
* __Phased Methylation__: Nanopolish + PRINCESS-subtool
* __QC Statistics__ for each step

---

## Installation
Princess was tested on CentOS release 6.7, Conda version 4.7.12 is installed:
for more information about [Installing Conda press here](https://bioconda.github.io/user/install.html#install-conda, "Install Conda")
To download same Conda version [here](https://repo.continuum.io/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh "Conda 4.7.12")*

1. After conda is installed. Snakemake should be installed and yaml
~~~
conda install snakemake=5.7.1
conda install pyyaml
~~~
2. Downloading PRINCESS  
~~~
git clone https://github.com/MeHelmy/princess.git
~~~

---

## Tutorial

To have an overview about princess write command `princess -h`.
You will have the following list of commands that we can use in princess.

~~~
usage: princess [-h] {all,align,sv,snv,variant,phase,overview} ...

Princess A framework for long-reads analysis.

optional arguments:
  -h, --help            show this help message and exit

Sub-commands:
  Valid sub-commands

  {all,align,sv,snv,variant,phase,overview}
    all                 This command will run the following: Align the reads.
                        Identify SVs Identify SNVs Phase both SNVs and SVs
    align               This command will use the input sequence files and
                        align them against the reference using either Minimap2
                        or NGMLR use -a to choose aligner otherwise Minimap2
                        will be used by default.
    sv                  This command will use bam file to identify SV using
                        Sniffles.
    snv                 This command will use bam file to identify SNVs usin
                        Clair.
    variant             This command will use bam file to identify SVs and
                        SNVs.
    phase               This command will use use reads to identify SNVs by
                        Clair and Phase them.
    overview            This command will show what steps will run.

princess version 0.01. use command -h for info.
~~~


Assume that we want only to run `snv` command, to know more about its option:

`princess snv -h`


~~~
usage: princess snv [-h] [-v] -d Working directory -r {ont,clr,ccs} [-l] [-u]
                    [-e] [-a {minimap,ngmlr}]
                    [-s sampleFiles [sampleFiles ...]] -f REF [-j JOBS]
                    [-g LOG_FILE] [-c CHRS [CHRS ...]] [-t]

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -d Working directory, --directory Working directory
                        Working directory.
  -r {ont,clr,ccs}, --ReadType {ont,clr,ccs}
                        Read techonlogy (Note: clr is not supported anymore by clair3)
  -l, --removeFiles     remove princess source script after running default:
                        False)
  -u, --UseConda        Use conda for running default: True)
  -e, --Cluster         Use cluster while runing default: True)
  -a {minimap,ngmlr}, --Aligner {minimap,ngmlr}
                        In case if you want to choose specific aligner
                        otherwise default will be used default: minimap)
  -s sampleFiles [sampleFiles ...], --samplesFiles sampleFiles [sampleFiles ...]
                        list of fatsa, fastq, or gz files.
  -f REF, --ref REF     The reference file will be used to align reads to.
  -j JOBS, --jobs JOBS  Number of running jobs default: 200 )
  -g LOG_FILE, --log LOG_FILE
                        Log file: PrincessLog.txt )
  -c CHRS [CHRS ...], --chr CHRS [CHRS ...]
                        Chromosomes list, if not specified we will use all
                        Chromosomes.
  -t, --filter          Filter identified SNVs using Princess algorithm
                        default: True)
~~~


~~~
princess all  -d ./princess_all -r ont -s reads.split00.fastq.gz reads.split01.fastq.gz  -f hs37d5_mainchr.fa
~~~

`-r` defines the reads type.  
`-s` samples that we would like to analyze.  
`-f` **full path** to the reference.  

*__Note__*  
I am assuming that the reference file is indexed, if not please use the following command.  
`samtools faidx hs37d5_mainchr.fa` as a result you will have `hs37d5_mainchr.fa.fai`.

Done!!

## Output

Princess will create these directories:
- align   contains directory [minimap or ngmlr] based on the aligner that was specified.
- sv      contains the structural variant file sv/minimap/sniffles.vcf
- snp     contains single nucleotide variant calls per chromosomes
- phased  contains phased variant
- stat    contains Statistics
- meth    contains methylation info (if user choose to run methylation) 

---

## Converting from PBS to Slurm
1- Please ensure that you modify the `cluster/cluster_config.yaml` to specify the appropriate long-running node. For example, you can set the long queue system as follows:
    `long: &long_queue long_queue` 
    Where long_queue is the queue system that can run for a long time. Similarly, you can set the short queue in the following way:
    `short: &short_queue short_queue`, . Please refer to your cluster system administrator for more details.  
2- Please, ensure that you changed `cluster/config.yaml` from `cluster-status: "pbs_status.py"` to `cluster-status: "slurm_status.py"`  
3- In the `cluster/key_mapping.yaml` file. Please, change `system: "pbs"` to `system: "slurm"`  
4- Finally, in the `cluster/cluster_config.yaml` file, I set CPU and memory to each job to suit my cluster.  
E.g.
```
minimap2:
  queue: *long_queue
  time: "72:00:00"
  nCPUs: "12"
  mem: 20G
```
Here, I am using 12 CPUs, 20G memory, and the job running time is "72:00:00" maximum (three days.). You may need to use different configuration based on the availability in you cluster. Please, refer to your system administrator for more details.

