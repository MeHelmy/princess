# import Lib
############
import os


############################



# Include all snakefiles sub-moduels
###################################
prefixed = [filename for filename in os.listdir('./modules') if filename.endswith(".snakefile")]
for f in prefixed:
    include: f
