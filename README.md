# Princess
![princess](./pictures/leia.jpg)

## What it can do

* __Mapping__:  NGMLR or Minimap2
* __SNVs__: Clair (successor of Clairvoyante)
* __SVs__: Sniffles
* __Phasing SNVs__: WhatsHap
* __Phasing SVs__: PRINCESS-subtool
* __Extend Phasing__: PRINCESS-subtool
* __Phased Methylation__: Nanopolish + PRINCESS-subtool
* __QC Statistics__ for each step

## Installation
*Here I'm assuming that conda is installed.  
Princess was tested on CentOS release 6.7, Conda version 4.7.12
for more information about [Installing Conda press here](https://bioconda.github.io/user/install.html#install-conda, "Install Conda")
To download same Conda version [here](https://repo.continuum.io/miniconda/Miniconda3-4.7.12-Linux-x86_64.sh "Conda 4.7.12")*

1. After conda is installed. Snakemake should be installed
~~~
conda install snakemake=5.7.1
~~~
2. Downloading PRINCESS  
~~~
git clone git@github.com:MeHelmy/princess.git
~~~
3. Install Clair, Training models, pypy, and intervaltree
~~~
cd princess
chmod +x install.sh
./install.sh
~~~
### __Done__

## Tutorial
~~~
cd princess
~~~
Then copy your data or soft link it.  
Example: assume we are working on two files in `/home/user/samples`
and this directory contains two samples `sampl1.fasta` and `sample2.fatsa`
~~~
To copy:
cp /home/user/samples/sample1.fasta ./
cp /home/user/samples/sample12.fasta ./
Or you can soft link them like:
ln -s /home/user/samples/sample1.fasta ./
ln -s /home/user/samples/sample12.fasta ./
~~~
As you can see the files extensions in the previous step was `fasta` Now we will update the field  
`sample_extension` in the `config.yaml` instead of `gz` which is the default to be `fasta` the end result will be: `sample_extension: "fasta"`  
same step if it was `fastq, fa, fq etc...`

Time to run princess


## Output
