# import Lib
############
import os, glob, ntpath


############################


# Config File
#############
configfile: "./config.yaml"
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
