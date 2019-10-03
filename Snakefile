# import Lib
############
import os, glob, ntpath


############################


# Config File
#############
if len(config) == 0:
  if os.path.isfile("./config.yaml"):
    configfile: "./config.yaml"
  else:
    sys.exit("Looks like there is no config.yaml file in " + os.getcwd() + " make sure there is one or at least specify one with the --configfile commandline parameter.")
#############



# Listing samples
#################

# GET SAMPLES EXTENSION
sample_extension = config['sample_extension'] if config['sample_extension'] else "gz"

# GET WORKING DIRECTORY DEFAULT IS CURRENT DIRECTORY
data_dir = config["data_dir"] if config['data_dir'] else os.getcwd()

# GET SAMPLES LIST
sample_list = config['sample_list'] if config['sample_list'] else [ntpath.basename(sample) for sample in glob.glob( data_dir +"/*."+ sample_extension)]




# Preparing conda environements.
###############################
ALIGN="envs/align.yaml"
SV="envs/sv.yaml"
WHATSHAP="envs/whatshap.yaml"


# Importing scripts
###################
rawcoverage_script = config['read_raw_coverage']


# Include all snakefiles sub-moduels
###################################
prefixed = ["./modules/"+filename for filename in os.listdir('./modules') if filename.endswith(".smk")]
for f in prefixed:
    include: f
###################################

# RULES
#######
onstart:
    shell("cat pictures/start.txt")

rule all:
    input: expand("sv/{aligner}/sniffles.vcf", aligner=config['aligner'])

## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##
onsuccess:
	shell("cat pictures/success.txt")

onerror:
	shell("cat pictures/fail.txt")
