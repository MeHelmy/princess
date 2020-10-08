#!/usr/bin/env python3


import os
import sys
import warnings
import subprocess


jobid = sys.argv[1]

out= subprocess.run(['bjobs','-noheader',jobid],stdout=subprocess.PIPE).stdout.decode('utf-8')

state = out.strip().split()[2]


map_state={"PEND":'running',
           "RUN":'running',
           "PROV":"running",
           "WAIT":'running',
           "DONE":'success',
           "":'success'}

print(map_state.get(state,'failed'))
