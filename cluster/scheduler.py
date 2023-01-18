#!/usr/bin/env python3


import sys, os
import subprocess
from subprocess import Popen, PIPE
import yaml


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def run_cmd(cmd):
    eprint("Running from subprocess"+str(cmd))
    try:
        subprocess.run(cmd, check=True, universal_newlines=True)
    except subprocess.CalledProcessError as e:
        logger.error("Error in subprocess:\n{}".format(e.returncode))

def convert_time_to_seconds(run_time):
    # Number of : in input:
    colons = run_time.count(':')
    if colons == 3:
        # dd:hh:mm:ss
        d, h, m, s = run_time.split(':')
        return str(day2sec(int(d)) + hours2sec(int(h)) + minutes2sec(int(m)) + int(s))
    elif colons == 2:
        # hh:mm:ss
        h, m, s = run_time.split(':')
        return str( hours2sec(int(h)) + minutes2sec(int(m)) + int(s))
    elif colons == 1:
        # mm:ss
        m, s = run_time.split(':')
        return str(minutes2sec(int(m)) + int(s))
    else:
        return run_time

def day2sec(days):
    return days * 24 * 60 * 60

def hours2sec(hours):
    return hours * 60 * 60

def minutes2sec(minutes):
    return minutes * 60

# let snakemake read job_properties
from snakemake.utils import read_job_properties



jobscript = sys.argv[1]

job_properties = read_job_properties(jobscript)


#default paramters defined in cluster_spec (accessed via snakemake read_job_properties)
cluster_param= job_properties["cluster"]


if job_properties["type"]=='single':
    cluster_param['name'] = "snakejob.{}".format(job_properties['rule'])
elif job_properties["type"]=='group':
    cluster_param['name'] = job_properties['groupid']
else:
    raise NotImplementedError(f"Don't know what to do with job_properties['type']=={job_properties['type']}")


# don't overwrite default parameters if defined in rule (or config file)
if ('threads' in job_properties) and ('threads' not in cluster_param):
    cluster_param["threads"] = job_properties["threads"]
for res in ['time','mem']:
    if (res in job_properties["resources"]) and (res not in cluster_param):
        cluster_param[res] = job_properties["resources"][res]

## TODO: Comment this line || test while using normal time 01:00:00:00
# time in hours
if "time" in cluster_param:
    # cluster_param["time"]=int(cluster_param["time"])*60
    cluster_param["time"]=convert_time_to_seconds(cluster_param["time"])
    # cluster_param["time"]=cluster_param["time"]


# check which system you are on and load command command_options
key_mapping_file=os.path.join(os.path.dirname(__file__),"key_mapping.yaml")
command_options=yaml.load(open(key_mapping_file),
                          Loader=yaml.BaseLoader)
system= command_options['system']
command= command_options[system]['command']

key_mapping= command_options[system]['key_mapping']

# construct command:
for  key in key_mapping:
    if key in cluster_param:
        command+=" "
        command+=key_mapping[key].format(cluster_param[key])

command+=' {}'.format(jobscript)

eprint("submit command: "+command)

# run_cmd(command.split(' '))
p = Popen(command.split(' '), stdout=PIPE, stderr=PIPE)
output, error = p.communicate()
if p.returncode != 0:
    raise Exception("Job can't be submitted\n"+output.decode("utf-8")+error.decode("utf-8"))
else:
    res= output.decode("utf-8")

    if system=='lsf':
        import re
        match = re.search(r"Job <(\d+)> is submitted", res)
        jobid = match.group(1)

    elif system=='pbs':
        jobid= res.strip().split('.')[0]

    else:
        jobid= int(res.strip().split()[-1])

    print(jobid)
