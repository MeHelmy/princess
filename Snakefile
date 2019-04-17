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
samples = [ntpath.basename(sample) for sample in glob.glob(config['sample_directory']+"/*")]
extension = samples[0].rsplit(".", 1)[-1]
sample_list = [ i.rsplit(".", 1)[0] for i in samples ]



# Preparing conda environements.
###############################
ALIGN="envs/align.yaml"




# Include all snakefiles sub-moduels
###################################
prefixed = ["./modules/"+filename for filename in os.listdir('./modules') if filename.endswith(".snakefile")]
for f in prefixed:
    include: f
###################################

# RULES
#######

rule all:
    input: expand("align/{aligner}/data.bam", aligner=config['aligner'])
