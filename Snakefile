
# import Lib
############
import os, glob, ntpath, math
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



# Config reference and chromosomes list
#######################################
REFERENCES = config["reference"]
chr_list = config['chrs']
# ref = config["ref"]
# chr_list = {}
#
# if  config["chr_list"] and ref:
#     for reference in ref:
#         if reference in config["chr_list"] and config["chr_list"][reference]:
#             chr_list[reference] = config["chr_list"][reference]
#         else:
#             if os.path.isfile(REFERENCES[reference]+".fai"):
#                 chr_names = []
#                 with open(REFERENCES[reference]+".fai", 'r') as data_in:
#                     for line in data_in:
#                         chr_names.append(str(line.split()[0]))
#             else:
#                 print("Please make sure that {ref}.fai exists.\nOtherwise run:\nsamtools faidx {ref}".format(ref=REFERENCES[reference]))
#                 exit(1)
#             # f = Fasta(REFERENCES[reference])
#             # chr_list[reference] = [chr_name for chr_name in f.keys()]
#             chr_list[reference] = chr_names

# chromosomes List splited to chunks
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



# Metytlation variables
#######################
ont_sample_dir = config['fast5_dir']
#############



# Preparing conda environements.
###############################
PRINCESS_ENV=os.getcwd()+"/envs/princess_env.yaml"
# CLAIR_ENV=os.getcwd()+"/envs/test-clair.yaml"
CLAIR_ENV=os.getcwd()+"/envs/clair_env.yaml"
#############



# Importing scripts
###################
rawcoverage_script = config['read_raw_coverage']
updat_sv = config['updat_sv']
#############



# Include all snakefiles sub-moduels
###################################
prefixed = ["./modules/"+filename for filename in os.listdir('./modules') if filename.endswith(".smk")]
for f in prefixed:
    include: f
###################################



# Building output
##################
final_output = []

if not config['methylation']:
    pass
elif config['methylation'] and all(value  for value in ont_sample_dir.values()):
    final_output.append(data_dir + "/meth/"+ aligner + "/methylation_calls.tsv")
else:
    sys.exit("Every ONT sample should have corresponding fast5 directory, please correct fast5_dir files in config.yaml")

if config['update_snps'] and config['paternal_snps'] and config['maternal_snps']:
    final_output.extend([data_dir + "/stat.txt", *expand(data_dir + "/phased/{aligner}/data_updated.vcf", aligner=config['aligner']),\
    *expand(data_dir + "/sv/{aligner}/sniffles_hp_updated.vcf", aligner=config['aligner'])])
else:
    final_output.extend([data_dir + "/stat.txt", *expand(data_dir + "/sv/{aligner}/sv_snp.vcf.gz",  aligner=config['aligner'])])
###################################


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
	shell("mkdir -p {data_dir}/snake_log &&\
    find . -maxdepth 1 -name 'snakejob.*' -type f -print0 | xargs -0r mv -t {data_dir}/snake_log &&\
    cat pictures/success.txt")

onerror:
	shell("mkdir -p {data_dir}/snake_log &&\
    find . -maxdepth 1 -name 'snakejob.*' -type f -print0 | xargs -0r mv -t {data_dir}/snake_log &&\
    cat pictures/fail.txt")
