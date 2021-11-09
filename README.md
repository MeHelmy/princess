# Princess
[![GitHub](https://img.shields.io/github/license/MeHelmy/princess)](https://opensource.org/licenses/MIT) ![GitHub last commit](https://img.shields.io/github/last-commit/MeHelmy/princess) 
---
Princess is a fast and scalable framework to detect and report haplotype resolved Single Nucleotide Variants (SNV) and Structural Variations (SVs) at scale. It can leverage your cluster environment to speed up the detection which starts with one or many fasta or fastq files.
Cite the code: [![DOI](https://zenodo.org/badge/179986953.svg)](https://zenodo.org/badge/latestdoi/179986953)


![princess](./pictures/leia.jpg)

## Princess

* __Mapping__:  NGMLR or Minimap2
* __SNVs__: Clair (successor of Clairvoyante)
* __SVs__: Sniffles
* __Phasing SNVs__: WhatsHap
* __Phasing SVs__: PRINCESS-subtool
* __Extend Phasing__: PRINCESS-subtool
* __Phased Methylation__: Nanopolish + PRINCESS-subtool
* __QC Statistics__ for each step

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
3. Install Clair, Training models, pypy, and intervaltree
~~~
cd princess
chmod +x install.sh
./install.sh
~~~


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
                        Read techonlogy
  -l, --removeFiles     remove princess source script after running default:
                        False)
  -u, --UseConda        Use conda for running default: True)
  -e, --Cluster         Use cluster while runing default: True)
  -a {minimap,ngmlr}, --Aligner {minimap,ngmlr}
                        In case if you want to choose specific aligner
                        otherwise default will be used default: minimap)
  -s sampleFiles [sampleFiles ...], --sampleFiles sampleFiles [sampleFiles ...]
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

### For methylation calling.
Methylation calling is a part from the `all` option.

```
optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit
  -d Working directory, --directory Working directory
                        Working directory.
  -r {ont,clr,ccs}, --ReadType {ont,clr,ccs}
                        Read techonlogy
  -l, --removeFiles     remove princess source script after running default: False)
  -u, --UseConda        Use conda for running default: True)
  -e, --Cluster         Use cluster while runing default: True)
  -a {minimap,ngmlr}, --Aligner {minimap,ngmlr}
                        In case if you want to choose specific aligner otherwise default will be used default: minimap)
  -s sampleFiles [sampleFiles ...], --sampleFiles sampleFiles [sampleFiles ...]
                        list of fatsa, fastq, or gz files.
  -f REF, --ref REF     The reference file will be used to align reads to.
  -j JOBS, --jobs JOBS  Number of running jobs default: 200 )
  -g LOG_FILE, --log LOG_FILE
                        Log file: PrincessLog.txt )
  -c CHRS [CHRS ...], --chr CHRS [CHRS ...]
                        Chromosomes list, if not specified we will use all Chromosomes.
  -t, --filter          Filter identified SNVs using Princess algorithm default: True)
  -m, --methylation     Identify methylation, mutually inclusive with -md default: False)
  -md Fast5 Directory, --methylationDirectory Fast5 Directory
                        Fast5 directory will be used to identify methylation mutually inclusive with option -m default: False)
```
By choosing the flag __`--methylation`__, Princess will call the methylation on the input data (ONT data), this option is inclusive with the option __`--methylationDirectory`__ which requires the fasta5 directory.

## Test case

We uploaded a HiFi compressed data file from the publically available HG002 data set.
The complete data set (High-fidelity 15kb long-read dataset of HG002, Ashkenazim Son.) is available [Here](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/)

To download the test data run the following command:
```
wget https://bcm.box.com/shared/static/sdml5d7csxprgu3cl5cve0lgv5jnrrlv --output-document  HiFi.fastq.gz
```
After download is finished you shall have a HiFi fastq file called `HiFi.fastq.gz`, to run the analysis test run the following command:
```
Full/Path/To/princess all  --directory $PWD/analysis --ReadType ccs --ref Path/To/Reference/genome.fa --jobs 7 --sampleFiles $PWD/HiFi.fastq.gz  --latency-wait 200 -p
```
all:           The command to run full analysis for other options please run `princess -h`  
---directory:  The out put directory it could be any name, use the full path, in my case the output is  same place.  
--ReadType:    Read type, the supported read types are clr, ccs, and ont.  
--ref:         Path to the reference please use samtools faidx with refernce before running Princess.  
--jobs:        Number of running jobs on cluster.  
--sampleFiles: Sample fastq file we downloaded, it could be more than one either compressed or not.  
--latency-wait 200 -p:  These are additional Snakemake option to wait 200 seconds before collecting output.  






## Output

Princess will create these directories:
- align   contains directory [minimap or ngmlr] based on the aligner that was specified.
- sv      contains the structural variant file sv/minimap/sniffles.vcf
- snp     contains single nucleotide variant calls per chromosomes
- phased  contains phased variant
- stat    contains Statistics
- meth    contains methylation info (if user choose to run methylation)      
