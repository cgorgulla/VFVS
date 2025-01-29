#!/usr/bin/env python3

# Copyright (C) 2019 Christoph Gorgulla
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
#
# This file is part of VirtualFlow.
#
# VirtualFlow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# VirtualFlow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with VirtualFlow.  If not, see <https://www.gnu.org/licenses/>.

# ---------------------------------------------------------------------------
#
# Description: Submit jobs
#
# ---------------------------------------------------------------------------


import os
import json
import boto3
import botocore
import re
import argparse
import sys
import time
from botocore.config import Config
from pathlib import Path
import shutil
import jinja2
import subprocess


def parse_config(filename):
    with open(filename, "r") as read_file:
        config = json.load(read_file)

    return config


def run_bash(config, current_workunit, jobline):

    # Get the template
    try:
        with open(config['bash_template']) as f:
            bash_template = jinja2.Template(f.read())
    except IOError as error:
        print(f"Cannot open the bash_template ({config['bash_template']})")
        raise error

    jobline_str = str(jobline)

    # how many jobs are there that we need to submit?
    subjobs_count = len(current_workunit['subjobs'])

    # Where are we putting this file?
    batch_workunit_base=Path(config['sharedfs_workunit_path']) / jobline_str
    batch_workunit_base.mkdir(parents=True, exist_ok=True)
    batch_submit_file = batch_workunit_base / "run.bash"

    template_values = {
        "job_letter": config['job_letter'],
        "job_name": config['job_name'],
        "threads_to_use": config['threads_to_use'],
        "array_start": "0",
        "array_end": (subjobs_count - 1),
        "workunit_id": jobline_str,
        "job_storage_mode": config['job_storage_mode'],
        "job_tgz": current_workunit['download_path'],
        "batch_workunit_base": batch_workunit_base.resolve().as_posix()
    }
    render_output = bash_template.render(template_values)

    try:
        with open(batch_submit_file, "w") as f:
            f.write(render_output)
    except IOError as error:
        print(f"Cannot write the workunit bash file ({batch_submit_file})")
        raise error

    try:
        with open("../workflow/run.bash", "a") as bash_out:
            bash_out.write(f"bash {batch_submit_file}\n")
    except IOError as error:
        print(f"Cannot write the workflow bash file")
        raise error

    cmd = [
        "bash",
        batch_submit_file
    ]

    try:
        ret = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.TimeoutExpired as err:
        raise Exception("timeout on submission to bash")

    ret.wait()

    current_workunit['status'] = {
        'vf_job_status': 'SUBMITTED',
        'job_name': f"vfvs-{config['job_letter']}-{jobline_str}",
        'job_id': ret.pid
    }

    # Slow ourselves down a bit
    time.sleep(0.1)

def submit_slurm(config, client, current_workunit, jobline):


    # Get the template
    try:
        with open(config['slurm_template']) as f:
            slurm_template = jinja2.Template(f.read())
    except IOError as error:
        print(f"Cannot open the slurm_template ({config['slurm_template']})")
        raise error


    jobline_str = str(jobline)

    # how many jobs are there that we need to submit?
    subjobs_count = len(current_workunit['subjobs'])

    # Where are we putting this file?
    batch_workunit_base=Path(config['sharedfs_workunit_path']) / jobline_str
    batch_workunit_base.mkdir(parents=True, exist_ok=True)
    batch_submit_file = batch_workunit_base / "submit.slurm"


    #
    template_values = {
        "job_letter": config['job_letter'],
        "job_name": config['job_name'],
        "threads_to_use": config['threads_to_use'],
        "array_start": "0",
        "array_end": (subjobs_count - 1),
        "slurm_cpus": config['slurm_cpus'],
        "slurm_account": config['slurm_account'],
        "slurm_partition": config['slurm_partition'],
        "workunit_id": jobline_str,
        "job_storage_mode": config['job_storage_mode'],
        "slurm_array_job_throttle": config['slurm_array_job_throttle'],
        "job_tgz": current_workunit['download_path'],
        "batch_workunit_base": batch_workunit_base.resolve().as_posix()
    }
    render_output = slurm_template.render(template_values)

    try:
        with open(batch_submit_file, "w") as f:
            f.write(render_output)
    except IOError as error:
        print(f"Cannot write the workunit slurm file ({batch_submit_file})")
        raise error


    # Run the slurm submit command and capture the job ID
    # information

    cmd = [
        "sbatch",
        batch_submit_file
    ]

    try:
        ret = subprocess.run(cmd, capture_output=True,
                         text=True, timeout=int(config['slurm_job_submission_timeout']))
    except subprocess.TimeoutExpired as err:
        raise Exception("timeout on submission to sbatch")

    while ret.returncode is None:
        time.sleep(0.1)

    if ret.returncode == 0:
        match = re.search(
                r'^Submitted batch job (?P<value>[-0-9]+)', ret.stdout, flags=re.MULTILINE)
        if(match):
            matches = match.groupdict()
            job_id = int(matches['value'])
        else:
            raise Exception("sbatch returned, but cannot parse output")
    else:
        raise Exception("sbatch did not return successfully")

    current_workunit['status'] = {
        'vf_job_status': 'SUBMITTED',
        'job_name': f"vfvs-{config['job_letter']}-{jobline_str}", 
        'job_id': job_id
    }

    # Slow ourselves down a bit
    time.sleep(0.1)


def submit_aws_batch(config, client, current_workunit, jobline):

    jobline_str = str(jobline)

    # how many jobs are there that we need to submit?
    subjobs_count = len(current_workunit['subjobs'])

    # AWS Batch doesn't allow an array job of only 1 -- so if it's one
    # we will launch 2, but the second will exit quickly since it has
    # no work

    if(subjobs_count == 1):
        subjobs_count = 2

    # Path to the data files
    object_store_input_path = f"{config['object_store_job_prefix_full']}"

    # Which queue to submit to
    batch_queue_number = ((jobline - 1) % int(config['aws_batch_number_of_queues'])) + 1


    try:
        response = client.submit_job(
            jobName=f'vfvs-{config["job_letter"]}-{jobline}',
            timeout={
                'attemptDurationSeconds': int(config["aws_batch_subjob_timeout"])
            },
            jobQueue=f"{config['aws_batch_prefix']}-queue{batch_queue_number}",
            arrayProperties={
                'size': subjobs_count
            },
            jobDefinition=f"{config['aws_batch_jobdef']}",
            containerOverrides={
                'resourceRequirements': [
                    {
                        'type': 'VCPU',
                        'value': config['aws_batch_subjob_vcpus'],
                    },
                    {
                        'type': 'MEMORY',
                        'value': config['aws_batch_subjob_memory'],
                    },
                ],
                'environment': [
                
                    {
                        'name': 'VFVS_RUN_MODE',
                        'value': "awsbatch"
                    },
                    {
                        'name': 'VFVS_JOB_STORAGE_MODE',
                        'value': "s3"
                    },
                    {
                        'name': 'VFVS_VCPUS',
                        'value': config['threads_to_use']
                    },
                    {
                        'name': 'VFVS_TMP_PATH',
                        'value': config['tempdir_default']
                    },
                    {
                        'name': 'VFVS_RUN_SEQUENTIAL',
                        'value': "0"
                    },
                    {
                        'name': 'VFVS_WORKUNIT',
                        'value': jobline_str
                    },
                    {
                        'name': 'VFVS_CONFIG_JOB_OBJECT',
                        'value': current_workunit['s3_download_path']
                    },
                    {
                        'name': 'VFVS_CONFIG_JOB_BUCKET',
                        'value': config['object_store_job_bucket']
                    },
                    
                ]
            }
        )

        current_workunit['status'] = "SUBMITTED"
        current_workunit['job_id'] = response['jobId']

        for subjob_id, subjob in current_workunit['subjobs'].items():
            subjob['status'] = "SUBMITTED"



    except botocore.exceptions.ClientError as error: 
        print("invalid")
        raise error

    # Slow ourselves down a bit
    time.sleep(0.1)

def process(config, start, stop):

    aws_config = Config(
        region_name=config['aws_region']
    )
    client = boto3.client('batch', config=aws_config)

    status = {}

    submit_type=config['batchsystem']


    # If this is bash, let's init a file
    if submit_type == "bash":
        with open("../workflow/run.bash", "w") as bash_out:
            bash_out.write("#!/usr/bin/env bash\n")

    # load the status file that is keeping track of the data
    with open("../workflow/status.json", "r") as read_file:
        status = json.load(read_file)
        workunits = status['workunits']

    for jobline in range(start, stop + 1):

        jobline_str = str(jobline)
        if jobline_str not in workunits:
            print(f"Jobline {jobline_str} was not found")
        else:
            current_workunit = workunits[jobline_str]

            print(f"Submitting jobline {jobline_str}... ")

            # Now see if any of them have been submitted before
            if 'status' in current_workunit:
                print("jobs were already submitted for this....")
            else:
                if(submit_type == "awsbatch"):
                    submit_aws_batch(config, client, current_workunit, jobline)
                elif(submit_type == "slurm"):
                    submit_slurm(config, client, current_workunit, jobline)
                elif(submit_type == "bash"):
                    run_bash(config, current_workunit, jobline)
                else:
                    print(f"Unknown submit type {submit_type}")




    # Output all of the information about the workunits into JSON so we can easily grab this data in the future
    
    print("Writing the json status file out")

    with open("../workflow/status.json.tmp", "w") as json_out:
        json.dump(status, json_out, indent=4)

    Path("../workflow/status.json.tmp").rename("../workflow/status.json")

    print("Making copy for status.submission as backup")
    shutil.copyfile("../workflow/status.json", "../workflow/status.submission.json")

    print("Done")

    if submit_type == "bash":
        os.chmod('../workflow/run.bash', 0o755)
        print("Run ../workflow/run.bash to run these job steps")


def main():

    if len(sys.argv) != 3:
        print('You must supply exactly two arguments -- the start jobline and the end jobline.')
        print('Joblines start from 1')
        sys.exit()

    config = parse_config("../workflow/config.json")
    process(config, int(sys.argv[1]), int(sys.argv[2]))


if __name__ == '__main__':
    main()
