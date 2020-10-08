#!/usr/bin/env python3

import sys
import subprocess
import xml.etree.cElementTree as ET

jobid = sys.argv[1]

try:
    res = subprocess.run("qstat -f -x {}".format(jobid), check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)

    xmldoc = ET.ElementTree(ET.fromstring(res.stdout.decode())).getroot()
    job_state = xmldoc.findall('.//job_state')[0].text

    if job_state == "C":
        exit_status = xmldoc.findall('.//exit_status')[0].text
        if exit_status == '0':
            print("success")
        else:
            print("failed")
    else:
        print("running")

except (subprocess.CalledProcessError, IndexError, KeyboardInterrupt) as e:
    print("failed")
