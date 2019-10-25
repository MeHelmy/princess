# import Lib
############
import os, glob, ntpath
from pyfaidx import Fasta


############################


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
REFERENCES = config["references"]
ref = config["ref"]
chr_list = {}
if  config["chr_list"] and ref:
    for reference in ref:
        if reference in config["chr_list"] and config["chr_list"][reference]:
            chr_list[reference] = config["chr_list"][reference]
        else:
            f = Fasta(config["references"][reference])
            chr_list[reference] = [chr_name for chr_name in f.keys()]
else:
    for reference in ref:
        f = Fasta(config["references"][reference])
        chr_list[reference] = [chr_name for chr_name in f.keys()]
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
CLAIR_ENV=os.getcwd()+"/envs/clair.yml"
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
