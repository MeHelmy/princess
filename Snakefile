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

# Include all snakefiles sub-moduels
###################################
prefixed = ["./modules/"+filename for filename in os.listdir('./modules') if filename.endswith(".snakefile")]
for f in prefixed:
    include: f
###################################

# Listing samples
#################
full_sample = [ntpath.basename(sample) for sample in glob.glob(config['sample_directory']+"/*")]
print(full_sample)
