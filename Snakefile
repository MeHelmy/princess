
# import Lib
############
import os, glob, ntpath, math, shutil
from snakemake.utils import min_version

############################

# Snake Version
###############
min_version("5.7.1")



# Config File
#############
# # if len(config) == 0:
if os.path.isfile("config.yaml"):
    configfile: "config.yaml"
else:
    sys.exit("Looks like there is no config.yaml file in " + os.getcwd() + " make sure there is one or at least specify one with the --configfile commandline parameter.")
#############



# Listing samples
#################
# GET SAMPLES EXTENSION
sample_extension = config['sample_extension'] if config['sample_extension'] else "gz"

# GET WORKING DIRECTORY DEFAULT IS CURRENT DIRECTORY
data_dir =  config["sample_directory"] if config['sample_directory'] else os.getcwd()

# GET SAMPLES LIST
sample_list = config['sample_list']
if not isinstance(sample_list, list):
    sample_list = sample_list.split()
#############



# Output sample name
###################
SAMPLE_NAME = config['sample_name']
#############


# Clean after success
####################
source_dir = config['delete_files']
samples_names = config['delete_samples']
def clean(source_dir, data_dir, samples_names):
    file_list = os.listdir(source_dir)
    if samples_names:
        for sample in samples_names: os.remove(os.path.join(data_dir, os.path.basename(sample)))
    for sample in file_list:
        if os.path.isfile(sample):
            os.remove(sample)
        else:
            shutil.rmtree(sample)
#############



# Config reference and chromosomes list
#######################################
REFERENCES = config["reference"]
chr_list = config['chrs']

# chromosomes List split to chunks
split_size = config['chr_split'] if config['chr_split'] and (config['chr_split'] >= 1000000) else 1000000
ref_index_file = REFERENCES+".fai"
chr_range = {}
with open(ref_index_file, 'r') as data_in:
    for line in data_in:
        chr, length = line.split()[0:2]
        if chr in chr_list:
            # Identify number of splits
            chr_split = int(length) // split_size
            chr_split = chr_split if chr_split > 1 else 1
            # step_value = int(length)//chr_split if chr_split > 0 else int(length)
            step_value = int(length)//chr_split
            ranges = list(range(0, int(length), step_value))
            if len(ranges) == chr_split + 1:
                ranges[-1] = int(length)
            else:
                ranges.append(int(length))
            ranges[0] = 1
            chr_range[chr] = ranges
#############



# Declare aligner
#################
aligner = config["aligner"]
#############



# Methylation variables
#######################
# ont_sample_dir = config['fast5_dir']
#############



# Preparing conda environments.
###############################
PRINCESS_ENV=os.getcwd()+"/envs/princess_env.yaml"
SNIFFLES_ENV=os.getcwd()+"/envs/sniffles.yaml"
#CLAIR_ENV=os.getcwd()+"/envs/clair3.yaml"
CLAIR_ENV=os.getcwd()+"/envs/clair3_no_depend.yaml"
MINIMAP2_ENV=os.getcwd()+"/envs/minimap2.yaml"
WHATSHAP_ENV=os.getcwd()+"/envs/whatshap.yaml"
VARIANT_ENV=os.getcwd()+"/envs/variant_tools.yaml"
READ_STAT_ENV=os.getcwd()+"/envs/pythonRun.yaml"
#############



# Importing scripts
###################
rawcoverage_script = config['read_raw_coverage']
updat_sv = config['updat_sv']
#############



# Include all snakemake files sub-modules
########################################
prefixed = ["./modules/"+filename for filename in os.listdir('./modules') if filename.endswith(".smk")]
for f in prefixed:
    include: f
###################################



# Building output
##################
final_output = []

if config['sample_list']:
    if not config['methylation']:
        pass
    # elif config['methylation'] and all(value  for value in ont_sample_dir.values()):
    elif config['methylation'] and  config['fast5_dir']:
        final_output.append(data_dir + "/result" + "/methylation.{}_calls_hap.tsv".format(aligner)) # DONE
    else:
        sys.exit("Every ONT sample should have corresponding fast5 directory, please correct fast5_dir files in config.yaml or use -md option")

    if config['update_snps'] and config['paternal_snps'] and config['maternal_snps']:
        final_output.extend([data_dir + "/result/.allReadsparental.{aligner}.txt".format(aligner=aligner)])
    else:
        final_output.extend([data_dir + "/result/.all.Reads.{}.txt".format(aligner)])
else:
    if config['update_snps'] and config['paternal_snps'] and config['maternal_snps']:
        final_output.extend([data_dir + "/result/.allReadsparental.{aligner}.txt".format(aligner=aligner)])
    else:
        final_output.extend([data_dir + "/result/.all.noReads.{}.txt".format(aligner)])


##############


# RULES
#######
onstart:
    shell("cat pictures/start.txt")

rule all:
    input: final_output

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
    clean(source_dir, data_dir, samples_names)
    if os.path.exists(os.path.join(data_dir, ".snakemake")):
        import shutil
        shutil.rmtree(os.path.join(data_dir, ".snakemake"), ignore_errors=True)
    shell("mkdir -p {data_dir}/snake_log &&\
    find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t {data_dir}/snake_log {{}}  \;  &&\
    cat {source_dir}/pictures/success.txt")


onerror:
    shell("mkdir -p {data_dir}/snake_log &&\
    find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t {data_dir}/snake_log {{}}  \;  &&\
    cat {source_dir}/pictures/fail.txt")
