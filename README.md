# Princess (beta version)
Princess is a fast and scalable framework to detect and report haplotype resolved Single Nucleotide Variants (SNV) and Structural Variations (SVs) at scale. It can leverage your cluster environment to speed up the detection which starts with one or many fasta or fastq files.

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

1. After conda is installed. Snakemake should be installed
~~~
conda install snakemake=5.7.1
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
Example: assume we are working on two files in `/home/user/samples` and this directory contains two fasta files `sample1.fasta` and `sample2.fasta`  
also we want the output in the directory `/home/user/result`.
Assume the Princess is installed in directory `/home/user/tools/princess`
~~~
/home/user/tools/princess/princess -c all -d /home/user/result -s /home/user/samples/sampl1.fasta  /home/user/samples/sampl2.fasta -r ont

# -r could be ont, clr, or ccs
~~~

### **Requirements**
You will need to adopt the `config.yaml` to suit your run. Especially:
~~~
ref: ['37']   # set it to the reference you want to run against either 37 or 38
references:   
  '37': "/home/user/reference/RCh37.fa"    # location of the indexed reference ex: /home/user/reference/RCh38.fa in the same directory you can use samtools faidx <ref.fasta>
  '38': "/home/user/reference/RCh38.fa"
~~~
If you want to use your cluster please define how many jobs should run:
~~~
cluster_jobs: 200 # default 200 change it to what you see convenient
~~~
Define your sequencing technology PacBio or ONT:
~~~
read_type: "pacbio" # choose from ont or pacbio
~~~
If you want to work on specific chromosomes, declare them in this variable remember to use the chromosomes names as it is in your reference file so If they are named 1,2 ... use this values other wise use chr1 chr2 ..:
~~~
chr_list:
  '37': [] # If I want to work on chromosomes 1 , 2 and 3 I shall update this field to be: [1,2,3]
  '38': [] #[chr1,chr2]
~~~
## Output

Princess will create these directories:
- align   contains directory [minimap or ngmlr] based on the aligner that was specified.
- sv      contains the structural variant file sv/minimap/sniffles.vcf
- snp     contains single nucleotide variant calls per chromosomes
- phased  contains phased variant
- stat    contains Statistics
- meth    contains methylation info (if user choose to run methylation)      
