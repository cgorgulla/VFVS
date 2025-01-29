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
# Description: Main runner for the individual workunits/job lines
#
# Revision history:
# 2021-06-29  Original version
# 2021-08-02  Added additional handling for case where there is only 
#             a single subjob in a job
# 2022-04-20  Adding support for output into parquet format
#
# ---------------------------------------------------------------------------

import tempfile
import tarfile
import gzip
import os
import json
import re
import boto3
import multiprocessing
import subprocess
import botocore
import logging
import time
import shutil
import hashlib
import pandas as pd
from pathlib import Path
from botocore.config import Config
from statistics import mean
from multiprocessing import Process
from multiprocessing import Queue
from queue import Empty
import math
import sys
import uuid

from typing import Dict, Tuple

logger = logging.getLogger(__name__)

def read_config_line(line: str) -> Tuple[str, str]:
    key, sep, value = line.strip().partition("=")
    return key.strip(), value.strip()

def load_config(config_path: str) -> Dict:
    with open(config_path) as fd:
        result = dict(read_config_line(line) for line in fd)

    for item in result:
        if '#' in result[item]:
            result[item] = result[item].split('#')[0]

    return result

def format_ligand(ligand_path: str, file_format: str) -> str:
    """Converts a ligand file to a different file format using the Open Babel tool.

        Args:
            ligand_ (str): The path to the input ligand file.
            new_format (str): The desired output format for the ligand file.
    
        Returns:
            None
    
        Raises:
            Exception: If the input file does not exist, or if the Open Babel tool is not installed.
    
        Examples:
            To convert a ligand file from mol2 format to pdbqt format:
            >>> convert_ligand_format('./ligands/ligand1.mol2', 'pdbqt')
    """
    ligand_path_as_list = ligand_path.split('.')
    current_format = ligand_path_as_list[-1]
    if current_format != file_format: 
        logger.info(f'Converting ligand file format to {file_format} using obabel.')
        ligand_path_as_list[-1] = file_format
        result = ''.join(ligand_path_as_list)
        os.system('obabel {} -O {}'.format(ligand_path, result))
        return result
        
    return ligand_path

def downloader(download_queue, unpack_queue, summary_queue, tmp_dir):


    botoconfig = Config(
       retries = {
          'max_attempts': 25,
          'mode': 'standard'
       }
    )

    s3 = boto3.client('s3', config=botoconfig)


    while True:
        try:
            item = download_queue.get(timeout=20.5)
        except Empty:
            continue

        if item is None:
            break


        item['temp_dir'] = tempfile.mkdtemp(prefix=f"{tmp_dir}/")
        item['local_path'] = f"{item['temp_dir']}/tmp.{item['ext']}"

        # Move the data either from S3 or a shared filesystem

        if('s3_download_path' in item['collection']):
            remote_path = item['collection']['s3_download_path']
            job_bucket = item['collection']['s3_bucket']

            try:
                with open(item['local_path'], 'wb') as f:
                    s3.download_fileobj(job_bucket, remote_path, f)
            except botocore.exceptions.ClientError as error:
                    reason = f"Failed to download from S3 {job_bucket}/{remote_path} to {item['local_path']}: ({error})"
                    logging.error(reason)

                    # log this to know we skipped
                    # put on the logging queue
                    summary_item = {
                        'type': "download_failed",
                        'log': {
                            'base_collection_key': item['collection_key'],
                            'reason': reason,
                            'dockings': item['collection']['dockings']
                        }
                    }
                    summary_queue.put(summary_item)
                    continue

        elif('sharedfs_path' in item['collection']):
            shutil.copyfile(Path(item['collection']['sharedfs_path']), item['local_path'])


        # Move it to the next step if there's space

        while unpack_queue.qsize() > 35:
            time.sleep(0.2)

        unpack_queue.put(item)




def untar(unpack_queue, collection_queue):
    print('Unpacker: Running', flush=True)

    while True:
        try:
            item = unpack_queue.get(timeout=20.5)
        except Empty:
            #print('unpacker: gave up waiting...', flush=True)
            continue

        # check for stop
        if item is None:
            break

        item['ligands'] = unpack_item(item)

        # Next step is to process this
        while collection_queue.qsize() > 35:
            time.sleep(0.2)

        collection_queue.put(item)



def unpack_item(item):

    ligands = {}

    os.chdir(item['temp_dir'])
    try:
        tar = tarfile.open(item['local_path'])
        for member in tar.getmembers():
            if(not member.isdir()):
                _, ligand = member.name.split("/", 1)

                if(ligand == ".listing"):
                    continue

                ligand_name = ligand.split(".")[0]

                ligands[ligand_name] = {
                    'path':  os.path.join(item['temp_dir'], item['collection']['collection_number'], ligand),
                    'base_collection_key': item['collection_key'],
                    'collection_key': item['collection_key']
                }

        tar.extractall()
        tar.close()
    except Exception as err:
        logging.error(
            f"ERR: Cannot open {item['local_path']} type: {str(type(err))}, err: {str(err)}")
        return None

    # Check if we have specific instructions
    if(item['collection']['mode'] == "prescreen_mode"):
        prescreen_ligands = {}

        # We should read the sparse file to know which ligands we actually need to keep
        with open(os.path.join(item['temp_dir'], item['collection']['collection_number'], ".listing"), "r") as read_file:
            for index, line in enumerate(read_file):
                line = line.strip()

                screen_collection_key, screen_ligand_name, screen_index = line.split(",")

                if(int(screen_index) >= item['collection']['prescreen_ligands_per_tranche']):
                    continue

                prescreen_ligands[screen_ligand_name] = ligands[screen_ligand_name]
                prescreen_ligands[screen_ligand_name]['collection_key'] = screen_collection_key

        return prescreen_ligands
    elif(item['collection']['mode'] == "named"):
        select_ligands = {}
        for ligand_name in item['collection']['ligands']:
            select_ligands[ligand_name] = ligands[ligand_name]
        return select_ligands
    else:
        return ligands





def collection_process(ctx, collection_queue, docking_queue, summary_queue):

    while True:
        try:
            item = collection_queue.get(timeout=20.5)
        except Empty:
            #print('unpacker: gave up waiting...', flush=True)
            continue

        # check for stop
        if item is None:
            break

        expected_ligands = 0
        completions_per_ligand = 0


        # How many do we run per ligand?
        for scenario_key in ctx['main_config']['docking_scenarios']:
            scenario = ctx['main_config']['docking_scenarios'][scenario_key]
            for replica_index in range(scenario['replicas']):
                completions_per_ligand += 1


            # generate directory for output
            scenario_directory = Path(item['temp_dir']) / "output" / scenario_key
            scenario_directory.mkdir(parents=True, exist_ok=True)


        # Process every ligand after making sure it is valid

        for ligand_key in item['ligands']:
            ligand = item['ligands'][ligand_key]

            coords = {}
            skip_ligand = 0

            with open(ligand['path'], "r") as read_file:
                for index, line in enumerate(read_file):

                    if (int(ctx['main_config']['run_atom_check']) == 1):
                        match = re.search(r'(?P<letters>\s+(B|Si|Sn)\s+)', line)
                        if(match):
                            matches = match.groupdict()
                            logging.error(
                                f"Found {matches['letters']} in {ligand}. Skipping.")
                            skip_reason = f"failed(ligand_elements:{matches['letters']})"
                            skip_reason_json = f"ligand includes elements: {matches['letters']})"
                            skip_ligand = 1
                            break

                    match = re.search(r'^ATOM', line)
                    if(match):
                        parts = line.split()
                        coord_str = ":".join(parts[5:8])

                        if(coord_str in coords):
                            logging.error(
                                f"Found duplicate coordinates in {ligand}. Skipping.")
                            skip_reason = f"failed(ligand_coordinates)"
                            skip_reason_json = f"duplicate coordinates"
                            skip_ligand = 1
                            break
                        coords[coord_str] = 1

                    # Checking for dynamic tranche filtering
                    if (int(ctx['main_config']['dynamic_tranche_filtering']) == 1):

                        # Checking if the line contains the tranche
                        match = re.search('Tranche:', line)
                        if (match):
                            # Obtaining the tranche
                            parts = line.split()
                            tranche  = parts[2].strip()

                            # Compiling the user's regex
                            try:
                                compiled_regex = re.compile(ctx['main_config']['dynamic_tranche_filtering_regex'])
                            except re.error as e:
                                print(f"Invalid regular expression: {e}")
                                return None

                            # Search for the pattern in the test string
                            match = compiled_regex.search(tranche)

                            # Checking if the tranche is not part of the regex
                            if not match:
                                print(f" Skipping ligand {ligand} due to dynamic tranche filtering. Tranche of ligand: {tranche}, regex:{ctx['main_config']['dynamic_tranche_filtering_regex']}")
                                #logging.error(f"Ligand {ligand} not contained in part of the library specified by regular expression. Skipping.")
                                skip_reason = f"failed(dynamic_tranche_filtering)"
                                skip_reason_json = f"dynamic tranche filtering"
                                skip_ligand = 1
                                break


            if(skip_ligand == 0):

                # We can submit this for processing
                ligand_attrs = get_attrs(ctx['main_config']['ligand_library_format'], ligand['path'], ctx['main_config']['print_attrs_in_summary'])
                submit_ligand_for_docking(ctx, docking_queue, ligand_key, ligand['path'], ligand['collection_key'], ligand['base_collection_key'], ligand_attrs, item['temp_dir'])
                expected_ligands += completions_per_ligand

            else:
                # log this to know we skipped
                # put on the logging queue
                summary_item = {
                    'type': "skip",
                    'log': {
                        'base_collection_key': ligand['base_collection_key'],
                        'collection_key': ligand['collection_key'],
                        'ligand_key': ligand_key,
                        'reason': skip_reason_json
                    }
                }
                summary_queue.put(summary_item)


        # Let the summary queue know that it can delete the directory after len(ligands)
        # number of ligands have been processed


        summary_item = {
            'type': "delete",
            'temp_dir': item['temp_dir'],
            'base_collection_key': item['collection_key'],
            'expected_completions': expected_ligands
        }

        summary_queue.put(summary_item)




def get_attrs(ligand_format, ligand_path, attrs = ['smi']):

    valid_formats = ['pdbqt', 'mol2', 'pdb', 'sdf']

    attributes = {}
    attributes_found = 0
    for key in attrs:
        attributes[key] = "N/A"

    if ligand_format in valid_formats:
        with open(ligand_path, "r") as read_file:
            for line in read_file:
                line = line.strip()

                match = re.search(r"SMILES_current:\s*(?P<smi>.*)$", line)
                if(match):
                    attributes['smi'] = match.group('smi')
                    attributes_found += 1

                match = re.search(r"SMILES:\s*(?P<smi>.*)$", line)
                if(match):
                    attributes['smi'] = match.group('smi')
                    attributes_found += 1

                match = re.search(r"SMILES_orig:\s*(?P<smi_orig>.*)$", line)
                if(match and 'smi_orig' in attrs):
                    attributes['smi_orig'] = match.group('smi_orig')
                    attributes_found += 1

                match = re.search(r"\* Heavy atom count:\s*(?P<hacount>.*)$", line)
                if(match and 'heavy_atom_count' in attrs):
                    attributes['heavy_atom_count'] = match.group('hacount')
                    attributes_found += 1

                match = re.search(r"\* logp_obabel:\s*(?P<logp_obabel>.*)$", line)
                if (match and 'logp_obabel' in attrs):
                    attributes['logp_obabel'] = match.group('logp_obabel')
                    attributes_found += 1

                match = re.search(r"\* obenergy:\s*(?P<obenergy>.*)$", line)
                if (match and 'obenergy' in attrs):
                    attributes['obenergy'] = match.group('obenergy')
                    attributes_found += 1

                if(attributes_found >= len(attrs)):
                    break

    return attributes


def obabel_check_energy(item, input_file, max_energy):

    #step_timer_start = time.perf_counter()

    #with open(input_file, "r") as read_file:
    #    lines = read_file.readlines()

    try:
        ret = subprocess.run(['obenergy', input_file], capture_output=True, text=True, timeout=item['obenergy_timeout'])
    except subprocess.TimeoutExpired as err:
        item['log']['reason'] = 'obenergy timed out'
        logging.error(item['log']['reason'])
        logging.error(f"stdout:\n{ret.stdout}\nstderr:{ret.stderr}\n")

    output_lines = ret.stdout.splitlines()

    kcal_to_kJ = 4.184

    if(len(output_lines) > 0):
        # obabel v2
        match = re.search(r"^TOTAL\s+ENERGY\s+\=\s+(?P<energy>-?\d+\.?\d*)\s+(?P<energy_unit>(kcal|kJ))", output_lines[-1])
        if(match):
            energy_value = float(match.group('energy'))
            if(match.group('energy_unit') == "kcal"):
                energy_value *= kcal_to_kJ

            if(energy_value <= float(max_energy)):
                return 1

            if(energy_value > float(max_energy)):
                item['log']['reason'] = f"energy-check failed for {item['collection_key']} {item['ligand_key']} {item['scenario_key']} ({energy_value} > {max_energy})"
                logging.error(item['log']['reason'])
                return 0

    item['log']['reason'] = 'No energy-check output'
    logging.error(item['log']['reason'])
    return 0


def posebusters_check(item, input_file, receptor_file, timeout):

    try:
        ret = subprocess.run(['bust', input_file, '-p', receptor_file, '--top-n', '1', '--outfmt', 'csv'], capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired as err:
        item['log']['reason'] = f"bust timed out for {item['collection_key']} {item['ligand_key']} {item['scenario_key']}"
        logging.error(item['log']['reason'])
        return 0

    output_lines = ret.stdout.splitlines()

    if(len(output_lines) > 0):
        match = output_lines[-1]
        posebusters_values = match.split(',')[2:]

        check = True
        for idx, posebusters_value in enumerate(posebusters_values):
            if posebusters_value == "False":
                check = False
                break

        if(check):
            return 1
        else:
            item['log']['reason'] = f"One or more docking pose checks failed for {item['collection_key']} {item['ligand_key']} {item['scenario_key']}"
            logging.error(item['log']['reason'])
            return 0

    item['log']['reason'] = f"No pose check output for {item['collection_key']} {item['ligand_key']} {item['scenario_key']}"
    logging.error(item['log']['reason'])
    return 0


def submit_ligand_for_docking(ctx, docking_queue, ligand_name, ligand_path, collection_key, base_collection_key, ligand_attrs, temp_dir):

    for scenario_key in ctx['main_config']['docking_scenarios']:
        scenario = ctx['main_config']['docking_scenarios'][scenario_key]

        for replica_index in range(scenario['replicas']):

            ligand_directory_directory = Path(temp_dir) / "output" / scenario_key / ligand_name / str(replica_index)
            ligand_directory_directory.mkdir(parents=True, exist_ok=True)

            docking_item = {
                'ligand_key': ligand_name,
                'ligand_path': ligand_path,
                'scenario_key': scenario_key,
                'collection_key': collection_key,
                'base_collection_key': base_collection_key,
                'config_path': scenario['config'],
                'program': scenario['program'],
                'program_long': scenario['program_long'],
                'input_files_dir':  os.path.join(ctx['temp_dir'], "vf_input", "input-files"),
                'timeout': int(ctx['main_config']['program_timeout']),
                'tools_path': ctx['tools_path'],
                'threads_per_docking': int(ctx['main_config']['threads_per_docking']),
                'temp_dir': temp_dir,
                'attrs': ligand_attrs,
                'output_dir': str(ligand_directory_directory),
                'run_energy_check': int(ctx['main_config']['run_energy_check']),
                'energy_max': int(ctx['main_config']['energy_max']),
                'obenergy_timeout': int(ctx['main_config']['obenergy_timeout']),
                'run_pose_check': int(ctx['main_config']['run_pose_check']),
                'pose_check_timeout': int(ctx['main_config']['pose_check_timeout'])
            }

            docking_queue.put(docking_item)

    while docking_queue.qsize() > 100:
            time.sleep(0.2)


def read_config_line(line):
    key, sep, value = line.strip().partition("=")
    return key.strip(), value.strip()


def docking_process(ctx, docking_queue, summary_queue):

    items_queue = {}

    while True:
        try:
            item = docking_queue.get(timeout=2)
        except Empty:
            # We should go ahead and process what is in the queue,
            # even though we might not have everything
            for scenario_key in items_queue:
                scenario = ctx['main_config']['docking_scenarios'][scenario_key]
                docking_process_batch(summary_queue, scenario, items_queue[scenario_key], ctx['temp_dir'])
                items_queue[scenario_key] = []
            continue

        # check for stop
        if item is None:
            # Process what is left in our queue
            for scenario_key in items_queue:
                scenario = ctx['main_config']['docking_scenarios'][scenario_key]
                docking_process_batch(summary_queue, scenario, items_queue[scenario_key], ctx['temp_dir'])
                items_queue[scenario_key] = []
            break

        scenario_key = item['scenario_key']
        if scenario_key not in items_queue:
            items_queue[scenario_key] = []

        scenario = ctx['main_config']['docking_scenarios'][scenario_key]
        items_queue[scenario_key].append(item)

        if len(items_queue[scenario_key]) >= scenario['batchsizes']:
            docking_process_batch(summary_queue, scenario, items_queue[scenario_key], ctx['temp_dir'])
            items_queue[scenario_key] = []



def docking_process_setup_common(item, docking_type, temp_dir):
    item['start_time'] = time.perf_counter()

    # temporary directory that will be wiped after this docking is complete

    item['uuid'] = str(uuid.uuid4())
    item['tmp_run_dir'] = Path(temp_dir) / "run" / item['uuid']
    item['tmp_run_dir'].mkdir(parents=True, exist_ok=True)

    # Make a copy of the input files so we have paths that make sense

    item['tmp_run_dir_input'] = Path(item['tmp_run_dir']) / "input-files"
    shutil.copytree(item['input_files_dir'], item['tmp_run_dir_input'])

    if docking_type == "batch":
        item['output_dir'] = Path(temp_dir) / "logs" / item['uuid']
        item['output_dir'].mkdir(parents=True, exist_ok=True)

    item['log_path'] = f"{item['output_dir']}/stdout"


def docking_process_clean_common(item):
    shutil.rmtree(item['tmp_run_dir'])
    item['seconds'] = time.perf_counter() - item['start_time']

def docking_process_batch(summary_queue, scenario, items, temp_dir):

    if(len(items) == 0):
        return

    batched_item = {
        'items': items,
        'program': scenario['program'],
        'execution_type': DOCKING_PROGRAMS[scenario['program']]['ligands'],
        'scenario_key': scenario['key']
    }

    for item in batched_item['items']:
        item['output_path'] = f"{item['output_dir']}/output"
        item['status'] = "failed"
        item['log_path'] = f"{item['output_dir']}/stdout"
        item['log'] = {
            'base_collection_key': item['base_collection_key'],
            'collection_key': item['collection_key'],
            'ligand_key': item['ligand_key'],
            'reason': ""
        }

    if batched_item['execution_type'] == "single":
        for item in batched_item['items']:
            ret = None
            start = time.perf_counter()

            docking_process_setup_common(item, "single", temp_dir)

            print(f"processing {item['ligand_key']}")

            try:
                cmd = program_runstring_array(item)
            except RuntimeError as err:
                logging.error(f"Invalid cmd generation for {item['ligand_key']} (program: '{item['program']}')")
                raise(err)

            try:
                ret = subprocess.run(cmd, capture_output=True,
                         text=True, cwd=item['tmp_run_dir_input'], timeout=item['timeout'])
            except subprocess.TimeoutExpired as err:
                item['log']['reason'] = f"timeout on {item['ligand_key']}"
                logging.error(item['log']['reason'])


            if ret != None:
                if ret.returncode == 0:
                    process_docking_completion(item, ret)
                else:
                    item['log']['reason'] = f"Non zero return code for {item['collection_key']} {item['ligand_key']} {item['scenario_key']}"
                    logging.error(item['log']['reason'])
                    logging.error(f"stdout:\n{ret.stdout}\nstderr:{ret.stderr}\n")


                # Place output into files
                with open(item['log_path'], "w") as output_f:
                    output_f.write(f"STDOUT:\n{ret.stdout}\n")
                    output_f.write(f"STDERR:\n{ret.stderr}\n")

                item['seconds'] = time.perf_counter() - item['start_time']

                print(f"processing {item['ligand_key']} - done in {item['seconds']}")

            docking_process_clean_common(item)

    elif batched_item['execution_type'] == "batch":
        ret = None

        docking_process_setup_common(batched_item, "batch", temp_dir)

        print(f"processing batch of {len(batched_item['items'])} items")

        # Mark which docking these ligands were associated with
        for item in batched_item['items']:
            with open(f"{item['output_dir']}/dock_uuid", "w") as output_f:
                output_f.write(f"{item['uuid']}\n")

        try:
            cmd = program_runstring_array_batch(batched_item)
        except RuntimeError as err:
            logging.error(f"Invalid cmd generation for batched execution (program: '{batched_item['program']}')")
            raise(err)

        try:
            ret = subprocess.run(cmd, capture_output=True,
                     text=True, cwd=batched_item['tmp_run_dir_input'], timeout=batched_item['timeout'])
        except subprocess.TimeoutExpired as err:
            reason = "Batched execution timed out"
            for item in batched_item['items']:
                item['log']['reason'] = reason
            logging.error(reason)

        if ret != None:
            if ret.returncode == 0:
                process_docking_completion_batch(batched_item, ret)
            else:
                reason = f"Non zero return code for batched execution"
                for item in batched_item['items']:
                    item['log']['reason'] = reason

                logging.error(reason)
                logging.error(f"stdout:\n{ret.stdout}\nstderr:{ret.stderr}\n")


            # Place output into files
            with open(batched_item['log_path'], "w") as output_f:
                output_f.write(f"STDOUT:\n{ret.stdout}\n")
                output_f.write(f"STDERR:\n{ret.stderr}\n")


        print(f"processing - done in {batched_item['seconds']}")
        docking_process_clean_common(batched_item)

    else:
        logging.error(f"Invalid ligand processing model for program {scenario['program']}")
        raise RuntimeError(f"Invalid ligand processing model for program {scenario['program']}")


    while summary_queue.qsize() > 200:
        time.sleep(0.2)

    batched_item['type'] = "docking_complete"
    summary_queue.put(batched_item)



def move_batch_logs(item, scenario_directories):

    scenario_dest = Path(scenario_directories[item['scenario_key']]) / "batch_exec" / item['uuid']
    batch_output = Path(item['output_dir'])
    for dir_file in batch_output.iterdir():
        shutil.move(str(dir_file), f"{str(scenario_dest)}/")
    shutil.rmtree(item['output_dir'])

def check_for_completion_of_collection_key(collection_completions, collection_key, scenario_directories):

    current_completions = collection_completions[collection_key]['current_completions']
    expected_completions = collection_completions[collection_key]['expected_completions']

    if current_completions == expected_completions:
        for scenario_key, scenario_dest in scenario_directories.items():
            scenario_directory = Path(collection_completions[collection_key]['temp_dir']) / "output" / scenario_key
            for dir_file in scenario_directory.iterdir():
                try:
                    shutil.move(str(dir_file), f"{scenario_dest}/")
                except shutil.Error as err:
                    logging.error(f"Could not move {dir_file} to {scenario_dest}, skipping...")
                    continue

        shutil.rmtree(collection_completions[collection_key]['temp_dir'])
        collection_completions.pop(collection_key, None)




def summary_process(ctx, summary_queue, upload_queue, metadata):

    print("starting summary process")
    start_time =  time.perf_counter()

    overview_data = {
        'metadata': metadata,
        'total_dockings': 0,
        'dockings_status': {
            'success': 0,
            'failed': 0
        },
        'skipped_ligands': 0,
        'skipped_ligand_list': [],
        'failed_list': [],
        'failed_downloads': 0,
        'failed_downloads_log': [],
        'failed_downloads_dockings': 0
    }

    summary_data = {}
    dockings_processed = 0
    collection_completions = {}
    scenario_directory = {}

    for scenario_key in ctx['main_config']['docking_scenarios']:
        summary_data[scenario_key] = {}

        scenario_directory[scenario_key] = Path(ctx['temp_dir']) / "output" / scenario_key / ctx['subjob_id']
        scenario_directory[scenario_key].mkdir(parents=True, exist_ok=True)


    while True:
        try:
            item = summary_queue.get()
        except Empty:
            #print('unpacker: gave up waiting...', flush=True)
            continue

        # check for stop
        if item is None:

            # Calculate how much time we spent
            overview_data['sec']  = time.perf_counter() - start_time

            # We need to generate the general overview data
            # even if we didn't process anything

            generate_overview_file(ctx, overview_data, upload_queue)

            # clean up anything extra that is around
            if(dockings_processed > 0):
                generate_summary_file(ctx, summary_data, upload_queue, ctx['temp_dir'])
                generate_output_archives(ctx, upload_queue, scenario_directory)
            break


        if(item['type'] == "delete"):
            if item['base_collection_key'] in collection_completions:
                collection_completions[item['base_collection_key']]['expected_completions'] = item['expected_completions']
                collection_completions[item['base_collection_key']]['temp_dir'] = item['temp_dir']
            else:
                collection_completions[item['base_collection_key']] = {
                    'expected_completions':item['expected_completions'],
                    'current_completions': 0,
                    'temp_dir': item['temp_dir']
                }

            check_for_completion_of_collection_key(collection_completions, item['base_collection_key'], scenario_directory)

        elif(item['type'] == "download_failed"):
            overview_data['failed_downloads_log'].append(item['log'])
            overview_data['failed_downloads'] += 1
            overview_data['failed_downloads_dockings'] += item['log']['dockings']


        elif(item['type'] == "skip"):

            overview_data['skipped_ligands'] += 1
            overview_data['skipped_ligand_list'].append(item['log'])

        elif(item['type'] == "docking_complete"):
            dockings_processed += 1

            for single_item in item['items']:
                overview_data['total_dockings'] += 1
                overview_data['dockings_status'][single_item['status']] += 1

                if(single_item['status'] == "success"):
                    summary_key = f"{single_item['ligand_key']}"

                    if summary_key not in summary_data[single_item['scenario_key']]:
                        summary_data[single_item['scenario_key']][summary_key] = {
                            'ligand': single_item['ligand_key'],
                            'collection_key': single_item['collection_key'],
                            'scenario': single_item['scenario_key'],
                            'scores': [ single_item['score'] ],
                            'attrs': single_item['attrs']
                        }

                    else:
                        summary_data[single_item['scenario_key']][summary_key]['scores'].append(single_item['score'])
                else:
                    # Log the failure
                    overview_data['failed_list'].append(single_item['log'])


                # See if this was the last completion for this collection_key

                if single_item['base_collection_key'] in collection_completions:
                    collection_completions[single_item['base_collection_key']]['current_completions'] += 1
                    check_for_completion_of_collection_key(collection_completions, single_item['base_collection_key'], scenario_directory)
                else:
                    collection_completions[single_item['base_collection_key']] = {
                        'expected_completions': -1,
                        'current_completions': 1,
                        'temp_dir': ""
                    }

            # Copy data if this was a batch execution
            if(item['execution_type'] == "batch"):
                move_batch_logs(item, scenario_directory)

        else:
            logging.error(f"received invalid summary completion {item['type']}")
            raise




def generate_tarfile(dir, tarname):
    os.chdir(str(Path(dir).parents[0]))

    with tarfile.open(tarname, "x:gz") as tar:
        tar.add(os.path.basename(dir))

    return os.path.join(str(Path(dir).parents[0]), tarname)


def generate_output_path(ctx, scenario_key, content_type, extension):


    if scenario_key != None:
        if(ctx['main_config']['job_storage_mode'] == "s3"):
            return f"{ctx['main_config']['object_store_job_prefix']}/{ctx['main_config']['job_name']}/{scenario_key}/{content_type}/{ctx['workunit_id']}/{ctx['subjob_id']}.{extension}"
        else:
            outputfiles_dir = Path(ctx['main_config']['sharedfs_output_files_path']) / scenario_key / content_type / str(ctx['workunit_id'])
            return f"{ctx['main_config']['sharedfs_output_files_path']}/{scenario_key}/{content_type}/{ctx['workunit_id']}/{ctx['subjob_id']}.{extension}"
    else:
        if(ctx['main_config']['job_storage_mode'] == "s3"):
            return f"{ctx['main_config']['object_store_job_prefix']}/{ctx['main_config']['job_name']}/{content_type}/{ctx['workunit_id']}/{ctx['subjob_id']}.{extension}"
        else:
            outputfiles_dir = Path(ctx['main_config']['sharedfs_output_files_path']) / content_type / str(ctx['workunit_id'])
            return f"{ctx['main_config']['sharedfs_output_files_path']}/{content_type}/{ctx['workunit_id']}/{ctx['subjob_id']}.{extension}"


def generate_output_archives(ctx, upload_queue, scenario_directories):

    for scenario_key, scenario_dir in scenario_directories.items():

        upload_tmp_dir = tempfile.mkdtemp(prefix=ctx['temp_dir'])

        # tar the file
        tar_name = f"{upload_tmp_dir}/{ctx['subjob_id']}.tar.gz"
        tar_gz_path = generate_tarfile(scenario_dir, tar_name)

        uploader_item = {
            'storage_type': ctx['main_config']['job_storage_mode'],
            'remote_path' : generate_output_path(ctx, scenario_key, "logs", "tar.gz"),
            's3_bucket' : ctx['main_config']['object_store_job_bucket'],
            'local_path': tar_name,
            'temp_dir': upload_tmp_dir
        }

        upload_queue.put(uploader_item)


def generate_overview_file(ctx, overview_data, upload_queue):

    upload_tmp_dir = tempfile.mkdtemp(prefix=ctx['temp_dir'])
    overview_json_location = f"{upload_tmp_dir}/overview.json.gz"

    with gzip.open(overview_json_location, "wt") as json_out:
        json.dump(overview_data, json_out)

    uploader_item = {
        'storage_type': ctx['main_config']['job_storage_mode'],
        'remote_path' : generate_output_path(ctx, None, "summary", "json.gz"),
        's3_bucket': ctx['main_config']['object_store_job_bucket'],
        'local_path': overview_json_location,
        'temp_dir': upload_tmp_dir
    }

    upload_queue.put(uploader_item)

def generate_summary_file(ctx, summary_data, upload_queue, tmp_dir):

    for scenario_key in ctx['main_config']['docking_scenarios']:

        if(len(summary_data[scenario_key]) == 0):
            break

        csv_ordering = ['ligand', 'collection_key', 'scenario', 'score_average', 'score_min']
        max_scores = 0

        # Need to run all of the averages
        for summary_key, summary_value in summary_data[scenario_key].items():

            if(len(summary_value['scores']) > max_scores):
                max_scores = len(summary_value['scores'])

            for index, score in enumerate(summary_value['scores']):
                summary_value[f"score_{index}"] = score

            summary_value['score_average'] = mean(summary_value['scores'])
            summary_value['score_min'] = min(summary_value['scores'])

            summary_value.pop('scores', None)

            # Update the attrs
            for attr_name, attr_value in summary_value['attrs'].items():
                summary_value[f'attr_{attr_name}'] = attr_value
                if f'attr_{attr_name}' not in csv_ordering:
                    csv_ordering.append(f'attr_{attr_name}')

            summary_value.pop('attrs', None)


            # For each collection tranche.. .explode them out
            collection_name, collection_number = summary_value['collection_key'].split("_")
            for letter_index, letter in enumerate(collection_name):
                summary_value[f'tranche_{letter_index}'] = letter


        # Ouput for each scenario

        if 'parquet' in ctx['main_config']['summary_formats']:

            # Now we can generate a parquet file with all of the data
            df = pd.DataFrame.from_dict(summary_data[scenario_key], "index")

            upload_tmp_dir = tempfile.mkdtemp(prefix=tmp_dir)

            uploader_item = {
                'storage_type': ctx['main_config']['job_storage_mode'],
                'remote_path' : generate_output_path(ctx, scenario_key, "parquet", "parquet"),
                's3_bucket': ctx['main_config']['object_store_job_bucket'],
                'local_path': f"{upload_tmp_dir}/summary.parquet",
                'temp_dir': upload_tmp_dir
            }


            df.to_parquet(uploader_item['local_path'], compression='gzip')
            upload_queue.put(uploader_item)


        if 'csv.gz' in ctx['main_config']['summary_formats']:

            for index in range(max_scores):
                csv_ordering.append(f"score_{index}")

            upload_tmp_dir = tempfile.mkdtemp(prefix=tmp_dir)

            with gzip.open(f"{upload_tmp_dir}/summary.txt.gz", "wt") as summmary_fp:

                # print header
                summmary_fp.write(",".join(csv_ordering))
                summmary_fp.write("\n")

                for summary_key, summary_value in summary_data[scenario_key].items():
                    ordered_summary_list = list(map(lambda key: str(summary_value.get(key, "")), csv_ordering))
                    summmary_fp.write(",".join(ordered_summary_list))
                    summmary_fp.write("\n")


            uploader_item_csv = {
                'storage_type': ctx['main_config']['job_storage_mode'],
                'remote_path' : generate_output_path(ctx, scenario_key, "csv", "csv.gz"),
                's3_bucket': ctx['main_config']['object_store_job_bucket'],
                'local_path': f"{upload_tmp_dir}/summary.txt.gz",
                'temp_dir': upload_tmp_dir
            }
            upload_queue.put(uploader_item_csv)



def upload_process(ctx, upload_queue):

    print('Uploader: Running', flush=True)

    botoconfig = Config(
       retries = {
          'max_attempts': 25,
          'mode': 'standard'
       }
    )

    s3 = boto3.client('s3', config=botoconfig)

    while True:
        try:
            item = upload_queue.get(timeout=20.5)
        except Empty:
            #print('unpacker: gave up waiting...', flush=True)
            continue

        # check for completion
        if item is None:

            # clean up anything extra that is around
            break

        # Save off the data we need
        # Basically.. if s3 then use boto3

        if(item['storage_type'] == "s3"):
            try:
                print(f"Uploading to {item['remote_path']}")
                response = s3.upload_file(item['local_path'], item['s3_bucket'], item['remote_path'])
            except botocore.exceptions.ClientError as e:
                logging.error(e)
                raise

        else:
            # if sharedfs.. .then just do a copy

            parent_directory = Path(Path(item['remote_path']).parent)
            parent_directory.mkdir(parents=True, exist_ok=True)

            shutil.copyfile(item['local_path'], item['remote_path'])


        # Get rid of the temp directory
        shutil.rmtree(item['temp_dir'])






def process_config(ctx):

    # Create absolute directories based on the other parameters

    ctx['main_config']['collection_working_path'] = os.path.join(
        ctx['temp_dir'], "collections")
    ctx['main_config']['output_working_path'] = os.path.join(
        ctx['temp_dir'], "output-files")

    if('summary_formats' not in ctx['main_config']):
        ctx['main_config']['summary_formats'] = {
                'csv.gz': 1
        }


    ctx['main_config']['docking_scenarios'] = {}

    for index, scenario in enumerate(ctx['main_config']['docking_scenario_names']):
        program_long = ctx['main_config']['docking_scenario_methods'][index]
        program = program_long

        # Special handing for smina* and gwovina*
        match = re.search(r'^(?P<program>smina|gwovina)', program_long)
        if(match):
            matches = match.groupdict()
            program = matches['program']
            logging.debug(f"Found {program} in place of {program_long}")
        else:
            logging.debug(f"No special match for '{program_long}'")

        ctx['main_config']['docking_scenarios'][scenario] = {
            'key': scenario,
            'config': os.path.join(ctx['temp_dir'], "vf_input", "input-files",
                                   ctx['main_config']['docking_scenario_names'][index],
                                   "config.txt"
                                   ),
            'program': program,
            'program_long': program_long,
            'replicas': int(ctx['main_config']['docking_scenario_replicas'][index]),
            'batchsizes': int(ctx['main_config']['docking_scenario_batchsizes'][index]),
        }



def get_workunit_from_s3(ctx, workunit_id, subjob_id, job_bucket, job_object, download_dir):
    # Download from S3

    download_to_workunit_file = "/".join([download_dir, "vfvs_input.tar.gz"])

    try:
        with open(download_to_workunit_file, 'wb') as f:
            ctx['s3'].download_fileobj(job_bucket, job_object, f)
    except botocore.exceptions.ClientError as error:
        logging.error(
            f"Failed to download from S3 {job_bucket}/{job_object} to {download_to_workunit_file}, ({error})")
        sys.exit(1)

    os.chdir(download_dir)

    # Get the file with the specific workunit we need to work on
    try:
        tar = tarfile.open(download_to_workunit_file)
        tar.extractall()
        file = tar.extractfile(f"vf_input/config.json")

        all_config = json.load(file)
        if(subjob_id in all_config['subjobs']):
            ctx['subjob_config'] = all_config['subjobs'][subjob_id]
        else:
            logging.error(f"There is no subjob ID with ID:{subjob_id}")
            # AWS Batch requires that an array job have at least 2 elements,
            # sometimes we only need 1 though
            if(subjob_id == "1"):
                sys.exit(0)
            else:
                sys.exit(1)
                raise RuntimeError(f"There is no subjob ID with ID:{subjob_id}")

        tar.close()
    except Exception as err:
        logging.error(
            f"ERR: Cannot open {download_to_workunit_file}. type: {str(type(err))}, err: {str(err)}")
        sys.exit(1)
        return None


    ctx['main_config'] = all_config['config']



def get_workunit_from_sharedfs(ctx, workunit_id, subjob_id, job_tar, download_dir):
    # Download from sharedfs

    download_to_workunit_file = "/".join([download_dir, "vfvs_input.tar.gz"])

    shutil.copyfile(job_tar, download_to_workunit_file)

    os.chdir(download_dir)

    # Get the file with the specific workunit we need to work on
    try:
        tar = tarfile.open(download_to_workunit_file)
        tar.extractall()
        file = tar.extractfile(f"vf_input/config.json")

        all_config = json.load(file)
        if(subjob_id in all_config['subjobs']):
            ctx['subjob_config'] = all_config['subjobs'][subjob_id]
        else:
            logging.error(f"There is no subjob ID with ID:{subjob_id}")
            raise RuntimeError(f"There is no subjob ID with ID:{subjob_id}")

        tar.close()
    except Exception as err:
        logging.error(
            f"ERR: Cannot open {download_to_workunit_file}. type: {str(type(err))}, err: {str(err)}")
        return None


    ctx['main_config'] = all_config['config']



def process_docking_completion(item, ret):
    item['status'] = "failed"

    if item['program'] not in DOCKING_PROGRAMS:
        raise RuntimeError(f"No completion function for {item['program']}")
    elif 'end' not in DOCKING_PROGRAMS[item['program']]:
        raise RuntimeError(f"No completion function for {item['program']}")
    else:
        DOCKING_PROGRAMS[item['program']]['end'](item, ret)


def process_docking_completion_batch(batch_item, ret):
    batch_item['status'] = "failed"

    if batch_item['program'] not in DOCKING_PROGRAMS:
        raise RuntimeError(f"No completion function for {batch_item['program']}")
    elif 'end' not in DOCKING_PROGRAMS[batch_item['program']]:
        raise RuntimeError(f"No completion function for {batch_item['program']}")
    else:
        DOCKING_PROGRAMS[batch_item['program']]['end'](batch_item, ret)


# Generate the run command for a given program

def program_runstring_array(task):

    cmd = []

    if task['program'] not in DOCKING_PROGRAMS:
        raise RuntimeError(f"No start function for {task['program']}")
    elif 'start' not in DOCKING_PROGRAMS[task['program']]:
        raise RuntimeError(f"No start function for {task['program']}")
    else:
        cmd = DOCKING_PROGRAMS[task['program']]['start'](task)

    return cmd

def program_runstring_array_batch(batch_item):

    cmd = []

    if batch_item['program'] not in DOCKING_PROGRAMS:
        raise RuntimeError(f"No start function for {batch_item['program']}")
    elif 'start' not in DOCKING_PROGRAMS[batch_item['program']]:
        raise RuntimeError(f"No start function for {batch_item['program']}")
    else:
        cmd = DOCKING_PROGRAMS[batch_item['program']]['start'](batch_item)

    return cmd


####### Docking program configurations

## MpSDockZN

def docking_start_MpSDockZN(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    task['MpSDockZN_tmp_file'] = os.path.join(task['tmp_run_dir'], "MpSDockZN_run.sh")
    
    with open(task['MpSDockZN_tmp_file'], 'a+') as f: 
        f.writelines(['export Chimera={}\n'.format(config_['chimera_path'])])
        f.writelines(['export DOCK6={}\n'.format(config_['dock6_path'])])
        f.writelines(['charge=`$Chimera/bin/chimera --nogui --silent {} charges.py`\n'.format(task['ligand_path'])])
        f.writelines(['antechamber -i {} -fi mol2 -o ligand_input.mol2 -fo mol2 -at sybyl -c gas -rn LIG -nc $charge -pf y\n'.format(task['ligand_path'])])
        f.writelines(['$DOCK6/bin/showbox < {} \n'.format(config_['dock6_path'])])
        f.writelines(['$DOCK6/bin/grid -i {} \n'.format(config_['grid_in'])])
        f.writelines(['./executables/MpSDock -i {} \n'.format(config_['dock_in'])])
        
    os.system('chmod 0700 {}'.format(task['MpSDockZN_tmp_file'])) # Assign execution permisions on script
    cmd = ['./{}'.format(task['MpSDockZN_tmp_file'])]
        
    return cmd

def docking_finish_MpSDockZN(item, ret): 
    try: 
        score_path = os.path.join(item['tmp_run_dir'], "receptor_input_docked_result.list")
        score_all = []
        with open(score_path, 'r') as f: 
            lines = f.readlines()
        for item in lines: 
            A = item.split(' ')
            A = [x for x in A if x != '']
            try: score_1, score_2, score_3, score_4, score_5 = float(A[0]), float(A[1]), float(A[2]),float(A[3]), float(A[4])
            except: continue 
            final_score = score_1 + score_2 + score_3 + score_4 + score_5
            score_all.append(final_score)
        item['score'] = min(score_all)   
 
        pose_path = os.path.join(item['tmp_run_dir'], "receptor_input_docked_result.mol2")
        shutil.move(pose_path, item['output_dir'])  
        item['status'] = "success"

    except: 
        logging.error("failed parsing")

## SEED
def docking_start_SEED(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    task['seed_tmp_file'] = os.path.join(task['tmp_run_dir'], "seed_run.sh")
    
    with open(task['seed_tmp_file'], 'a+') as f: 
        f.writelines(["charge=`{}/bin/chimera --nogui --silent {} charges.py`\n".format(config_['chimera_path'], task['ligand_path'])])
        f.writelines(["antechamber -i {} -fi mol2 -o ligand_gaff.mol2 -fo mol2 -at gaff2 -c gas -rn LIG -nc $charge -pf y\n".format(task['ligand_path'])])
        f.writelines(["python {} ligand_gaff.mol2 ligand_gaff.mol2 ligand_seed.mol2\n".format(task['mol2seed4_receptor_script'])])
        f.writelines(["{}/bin/chimera --nogui {} dockprep.py \n".format(config_['chimera_path'], config_['receptor'])])
        f.writelines(["python {} receptor.mol2 receptor.mol2 receptor_seed.mol2\n".format(task['mol2seed4_receptor_script'])])
        f.writelines(["\t\t-out:suffix out\n"])
        f.writelines(['{}/seed4 {} > log'.format(task['tools_path'], config_['seed_inp_file'])])
        
    os.system('chmod 0700 {}'.format(task['seed_tmp_file'])) # Assign execution permisions on script
    cmd = ['./{}'.format(task['seed_tmp_file'])]
        
    return cmd

def docking_finish_SEED(item, ret): 
    try: 
        score_path = os.path.join(item['tmp_run_dir'], "seed_best.dat")
        with open(score_path, 'r') as f: 
            lines = f.readlines()
        docking_score = float([x for x in lines[1].split(' ') if x != ''][4])
        item['score'] = min(docking_score)   
 
        pose_path = os.path.join(item['tmp_run_dir'], "ligand_seed_best.mol2")
        shutil.move(pose_path, item['output_dir'])  
        item['status'] = "success"

    except: 
        logging.error("failed parsing")


## HDock
def docking_start_HDock(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    task['HDock_tmp_file'] = os.path.join(task['tmp_run_dir'], "HDock_run.sh")
    
    with open(task['HDock_tmp_file'], 'a+') as f: 
        f.writelines(['{}/hdock {} {} -out Hdock.out'.format(task['tools_path'], config_['receptor'], task['lig_path'])])
        f.writelines(['{}/createpl Hdock.out top100.pdb -nmax 1 -complex -models'.format(task['tools_path']) ])
        
    os.system('chmod 0700 {}'.format(task['HDock_tmp_file'])) # Assign execution permisions on script
    cmd = ['./{}'.format(task['HDock_tmp_file'])]
        
    return cmd

def docking_finish_HDock(item, ret): 
    try: 
        pose_path = os.path.join(item['tmp_run_dir'], "model_1.pdb")        
        with open(pose_path, 'r') as f: 
            lines = f.readlines()
        docking_score = float(lines[3].split(' ')[-1])
        
        shutil.move(pose_path, item['output_dir'])  
        item['status'] = "success"
        item['score'] = docking_score   
    except: 
        logging.error("failed parsing")


## dock6
def docking_start_dock6(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Generate dock6 input file: 
    with open(os.path.join(task['tmp_run_dir'], "dock.in"), 'w') as f: 
        f.writelines(['conformer_search_type                                        flex\n'])
        f.writelines(['user_specified_anchor                                        no\n'])
        f.writelines(['limit_max_anchors                                            no\n'])
        f.writelines(['min_anchor_size                                              40\n'])
        f.writelines(['pruning_use_clustering                                       yes\n'])
        f.writelines(['pruning_max_orients                                          100\n'])
        f.writelines(['pruning_clustering_cutoff                                    100\n'])
        f.writelines(['pruning_conformer_score_cutoff                               25.0\n'])
        f.writelines(['pruning_conformer_score_scaling_factor                       1.0\n'])
        f.writelines(['use_clash_overlap                                            no\n'])
        f.writelines(['write_growth_tree                                            no\n'])
        f.writelines(['use_internal_energy                                          yes\n'])
        f.writelines(['internal_energy_cutoff                                       100.0\n'])
        f.writelines(['ligand_atom_file                                             {}\n'.format(task['ligand_path'])])
        f.writelines(['limit_max_ligands                                            no\n'])
        f.writelines(['receptor_site_file                                           {}\n'.format(config_['receptor_site_file'])])
        f.writelines(['max_orientations                                             500\n'])
        f.writelines(['chemical_matching                                            no\n'])
        f.writelines(['use_ligand_spheres                                           no\n'])
        f.writelines(['bump_filter                                                  no\n'])
        f.writelines(['score_molecules                                              yes\n'])
        f.writelines(['contact_score_primary                                        no\n'])
        f.writelines(['contact_score_secondary                                      no\n'])
        f.writelines(['grid_score_primary                                           yes\n'])
        f.writelines(['grid_score_secondary                                         no\n'])
        f.writelines(['grid_score_rep_rad_scale                                     1\n'])
        f.writelines(['grid_score_vdw_scale                                         1\n'])
        f.writelines(['grid_score_grid_prefix                                       grid\n'])
        f.writelines(['dock3.5_score_secondary                                      no\n'])
        f.writelines(['continuous_score_secondary                                   no\n'])
        f.writelines(['footprint_similarity_score_secondary                         no\n'])
        f.writelines(['pharmacophore_score_secondary                                no\n'])
        f.writelines(['descriptor_score_secondary                                   no\n'])
        f.writelines(['gbsa_zou_score_secondary                                     no\n'])
        f.writelines(['gbsa_hawkins_score_secondary                                 no\n'])
        f.writelines(['SASA_score_secondary                                         no\n'])
        f.writelines(['amber_score_secondary                                        no\n'])
        f.writelines(['minimize_ligand                                              yes\n'])
        f.writelines(['minimize_anchor                                              yes\n'])
        f.writelines(['minimize_flexible_growth                                     yes\n'])
        f.writelines(['use_advanced_simplex_parameters                              no\n'])
        f.writelines(['simplex_max_cycles                                           1\n'])
        f.writelines(['simplex_score_converge                                       0.1\n'])
        f.writelines(['simplex_cycle_converge                                       1.0\n'])
        f.writelines(['simplex_trans_step                                           1.0\n'])
        f.writelines(['simplex_rot_step                                             0.1\n'])
        f.writelines(['simplex_tors_step                                            10.0\n'])
        f.writelines(['simplex_anchor_max_iterations                                500\n'])
        f.writelines(['simplex_grow_max_iterations                                  500\n'])
        f.writelines(['simplex_grow_tors_premin_iterations                          0\n'])
        f.writelines(['simplex_random_seed                                          0\n'])
        f.writelines(['simplex_restraint_min                                        no\n'])
        f.writelines(['atom_model                                                   all\n'])
        f.writelines(['vdw_defn_file                                                {}/parameters\n'.format(config_['dock6_path'])])
        f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(config_['dock6_path'])])
        f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(config_['dock6_path'])])
        f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(config_['dock6_path'])])
        f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(config_['dock6_path'])])
        f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
        f.writelines(['write_orientations                                           no\n'])
        f.writelines(['num_scored_conformers                                        1\n'])
        f.writelines(['rank_ligands                                                 no\n'])
    
    cmd = ['{}/bin/dock6'.format(config_['dock6_path']), '-i', os.path.join(task['tmp_run_dir'], "dock.in")]

    return cmd

def docking_finish_dock6(item, ret): 
    try: 
        dock_file = [x for x in os.listdir(item['tmp_run_dir']) if 'ligand_out' in x]
        dock_file = [x for x in dock_file if 'mol2' in x][0]
        os.system('cp {} {}'.format(dock_file, item['output_dir']))
        
        # Save the results: 
        with open('./ligand_out_scored.mol2', 'r') as f: 
            lines = f.readlines()
        docking_score = float(lines[2].split(' ')[-1])           

        item['score'] = docking_score   
        item['status'] = "success"

    except: 
        logging.error("failed parsing")


## rosetta-ligand
def docking_start_rosetta_ligand(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    task['rosetta_tmp_file'] = os.path.join(task['tmp_run_dir'], "rosetta_run.sh")
    
    with open(task['rosetta_tmp_file'], 'a+') as f: 
        f.writelines(['export ROSETTA={}\n'.format(config_['ROSETTA_location'])])
        f.writelines(['obabel {} -O conformers.sdf --conformer --score rmsd --writeconformers --nconf 30\n'.format(task['ligand_path'])])
        f.writelines(['$ROSETTA/source/scripts/python/public/molfile_to_params.py -n LIG -p LIG --conformers-in-one-file conformers.sdf\n'])
        f.writelines(['cat {} LIG.pdb > complex.pdb\n'.format(task['receptor'])])
        f.writelines(['echo "END" >> complex.pdb\n'])
        f.writelines(["$ROSETTA/source/bin/rosetta_scripts.default.linuxgccrelease  \\\n"])
        f.writelines(["	-database $ROSETTA/database \\\n"])
        f.writelines(["\t@ options \\\n"])
        f.writelines(["\t\t-parser:protocol {} \\\n".format(config_['dock_xml_file_loc'])])
        f.writelines(["\t\t-parser:script_vars X={} Y={} Z={} \\\n".format(config_['center_x'], config_['center_y'], config_['center_z'])])
        f.writelines(["\t\t-in:file:s complex.pdb \\\n"])
        f.writelines(["\t\t-in:file:extra_res_fa LIG.params \\\n"])
        f.writelines(["\t\t-out:nstruct 10 \\\n"])
        f.writelines(["\t\t-out:level {} \\\n".format(config_['exhaustiveness'])])
        f.writelines(["\t\t-out:suffix out\n"])
        

    os.system('chmod 0700 {}'.format(task['rosetta_tmp_file'])) # Assign execution permisions on script
    cmd = ['./{}'.format(task['rosetta_tmp_file'])]
        
    return cmd

def docking_finish_rosetta_ligand(item, ret): 
    try: 
        
        docking_score_path = os.path.join(item['tmp_run_dir'], "scoreout.sc")

        with open(docking_score_path, 'r') as f: 
            lines = f.readlines()
        lines = lines[2: ]
        docking_scores = []
        for item in lines: 
            A = item.split(' ')
            A = [x for x in A if x!='']
            docking_scores.append(float(A[44]))
        
        item['score'] = min(docking_scores)   
        
        out_files = [x for x in os.listdir(item['tmp_run_dir']) if 'complexout' in x][0] 
        docking_out_path = os.path.join(item['tmp_run_dir'], out_files)
        shutil.move(docking_out_path, item['output_dir'])             
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
    return 


## LigandFit
def docking_start_Molegro(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    mvdscript_loc = os.path.join(task['tmp_run_dir'], "docking.mvdscript")
    
    with open(mvdscript_loc, 'a+') as f: 
        f.writelines(['// Molegro Script Job.\n\n'])
        f.writelines(['IMPORT Proteins;Waters;Cofactors FROM {}\n\n'.format(config_['receptor'])])
        f.writelines(['PREPARE Bonds=IfMissing;BondOrders=IfMissing;Hydrogens=IfMissing;Charges=Always; TorsionTrees=Always\n\n'])
        f.writelines(['IMPORT All FROM ligands.mol2\n\n'])
        f.writelines(['SEARCHSPACE radius=12;center=Ligand[0]\n\n'])
        f.writelines(['DOCK Ligand[1]\n\n\n'])
        f.writelines(['EXIT'])
            
            
    task['molegro_tmp_file'] = os.path.join(task['tmp_run_dir'], "run_.sh")
    
    with open(task['molegro_tmp_file'], 'w') as f: 
        f.writelines('export Molegro={}\n'.format(config_['molegro_location']))
        f.writelines('cat {} {} > ligands.mol2\n'.format(config_['ref_ligand'], task['ligand_path']))
        f.writelines('$Molegro/bin/mvd docking.mvdscript -nogui\n')

    os.system('chmod 0700 {}'.format(task['molegro_tmp_file'])) # Assign execution permisions on script
    cmd = ['./{}'.format(task['molegro_tmp_file'])]
        
    return cmd

def docking_finish_Molegro(item, ret): 
    try: 
        
        cmd_run = ret.stdout.split('\n')[-2]
        cmd_run = [x for x in cmd_run if 'Pose:' in x]
        scores = []
        for item in cmd_run: 
            scores.append( float(item.split('Energy')[-1].split(' ')[1][:-2]) )
        item['score'] = min(scores)
        
        docking_out_file = os.path.join(item['tmp_run_dir'])
        docking_out_file = [x for x in os.listdir(docking_out_file) if 'mol2' in x][0]
        docking_out_file = os.path.join(item['tmp_run_dir'], docking_out_file)

        shutil.move(docking_out_file, item['output_dir'])             
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
    return 

## FitDock
def docking_start_FitDock(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
    
    cmd = [
            f"{task['tools_path']}/FitDock",
            '-Tprot', config_['receptor_template'],
            '-Tlig', config_['ligand_reference'],
            '-Qprot', config_['receptor'],
            '-Qlig', task['ligand_path'], 
            '-ot', item['tmp_run_dir']+'ot.mol2', 
            '-os', item['tmp_run_dir']+'os.mol2', 
            '-o', item['tmp_run_dir']+'o.mol2'
          ]
        
    return cmd

def docking_finish_FitDock(item, ret): 
    try: 
        
        docking_score_file = os.path.join(item['tmp_run_dir'], "out.log")
        docking_out_file = os.path.join(item['tmp_run_dir'], "o.mol2")
                
        with open(docking_score_file, 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if 'Binding Score after  EM' in x]
        docking_score = float(lines[0].split(' ')[-2])
        item['score'] = min(docking_score)

        shutil.move(docking_out_file, item['output_dir'])             
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
    return 


## Flexx
def docking_start_flexx(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    output_file_path = item['output_dir'] + task['lig_path'].split('/')[-1].replace('.mol2', '.sdf')
    
    cmd = [f"{task['tools_path']}/flexx", 
           '-i', task['lig_path'], 
           '-o', output_file_path, 
           '-p', config_['receptor'], 
           '-r', config_['ref_ligand']
          ]
    
    return cmd

def docking_finish_flexx(item, ret): 
    try: 
        output_file_path = item['output_dir'] + item['lig_path'].split('/')[-1].replace('.mol2', '.sdf')
        
        with open(output_file_path, 'r') as f: 
            lines = f.readlines()    
        
        for i,item in enumerate(lines): 
            docking_scores = []
            if '>  <docking-score>' in item : 
                docking_score = float(lines[i+1])
                docking_scores.append(docking_score)
        
        item['score'] = min(docking_scores)
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
    return 


## LightDock
def docking_start_LightDock(task): 

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
    
    task['lightdock_tmp_file'] = os.path.join(task['tmp_run_dir'], "run_.sh")
    
    with open(task['lightdock_tmp_file'], 'w') as f: 
        f.writelines('{}/bin/lightdock3_setup.py {} {} --noxt --noh --now -anm\n'.format(config_['lightdock_path'], config_['receptor'], task['ligand_path']))
        f.writelines('{}/bin/lightdock3.py setup.json 100 -c 1 -l 0\n')
        f.writelines('{}/bin/lgd_generate_conformations.py {} {} swarm_0/gso_100.out {}\n'.format(config_['lightdock_path'], config_['receptor'], task['ligand_path'], config_['exhaustiveness']))

    os.system('chmod 0700 {}'.format(task['lightdock_tmp_file'])) # Assign execution permisions on script
    
    cmd = ['./{}'.format(task['lightdock_tmp_file'])]
    
    return cmd

def docking_finish_LightDock(item, ret): 
    try: 
        
        docking_score_file = os.path.join(item['tmp_run_dir'], "swarm_0", "gso_100.out")
        
        with open(docking_score_file, 'r') as f: 
            lines = f.readlines()
        lines = lines[1: ]
        scoring = []
        for item in lines: 
            A = item.split(' ')
            scoring.append(float(A[-1]))

        docking_pose_file = os.path.join(item['tmp_run_dir'], "swarm_0")
        complex_file = [x for x in os.listdir(docking_pose_file) if '.pdb' in x][0]
        docking_pose_file = os.path.join(item['tmp_run_dir'], "swarm_0", complex_file)
        shutil.move(docking_pose_file, item['output_dir'])     
        
        item['score'] = min(scoring)
        item['status'] = "success"
        
    except: 
        logging.error("failed parsing")
    return 

## RLDock
def docking_start_rldock(task): 
    cpus_per_program = str(task['threads_per_docking'])

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = [
            f"{task['tools_path']}/RLDOCK",
            '--i', config_['receptor'],
            '--l', task['ligand_path'],
            '-c', config_['exhaustiveness'],
            '-n', cpus_per_program, 
            '-s', config_['spheres_file_path']
        ]
    return cmd

def docking_finish_rldock(item, ret): 
    try: 
        docking_pose = os.path.join(item['tmp_run_dir_input'], "output_cluster.mol2")

        with open(docking_pose, 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if '# Total_Energy:' in x]
        docking_scores = []
        for item in lines: 
            docking_scores.append(float(item.split(' ')[-1]))

        shutil.move(docking_pose, item['output_dir'])        
        item['score'] = min(docking_scores)
        item['status'] = "success"

    except: 
        logging.error("failed parsing")
    
    return 

## Autodock koto
def docking_start_autodock_koto(task): 
    cpus_per_program = str(task['threads_per_docking'])

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = [
            f"{task['tools_path']}/AutoDock-Koto",
            '--receptor', config_['receptor'],
            '--ligand', task['ligand_path'],
            '--cpu', cpus_per_program,
            '--exhaustiveness', config_['exhaustiveness'],
            '--center_x', '{}'.format(config_['center_x']),
            '--center_y', '{}'.format(config_['center_y']),
            '--center_z', '{}'.format(config_['center_z']),
            '--size_x',   '{}'.format(config_['size_x']),
            '--size_y',   '{}'.format(config_['size_y']),
            '--size_z',   '{}'.format(config_['size_z']),
            '--out', task['output_path']
        ]
    return cmd

def docking_finish_autodock_koto(item, ret): 
    try: 
        docking_out = ret.stdout
        A = docking_out.split('\n')
        docking_score = []
        for item in A: 
            line_split = item.split(' ')
            line_split = [x for x in line_split if x != '']
            if len(line_split) == 4: 
                try: 
                    vr_1 = float(line_split[0])
                    vr_2 = float(line_split[1])
                    vr_3 = float(line_split[2])
                    vr_4 = float(line_split[3])
                    docking_score.append(vr_2)
                except: continue
            item['score'] = min(docking_score)
            item['status'] = "success"        
    except: 
        logging.error("failed parsing")
    
    return 


# CovDock
def docking_start_covdock(task):
    script_path = os.path.join(task['tmp_run_dir'], 'covdock_script.py')
    with open(script_path, 'w') as script_file:
        script_file.write(f"""\
        from schrodinger import structure
        from schrodinger.job import jobcontrol
        from schrodinger.application.covdock import covdock
        from schrodinger.application.ligprep import LigprepJob, LigprepSettings
        
        config_ = {task['config']}
        smi = '{task['ligand_path']}'
        receptor = config_['receptor']
        center_x = config_['center_x']
        center_y = config_['center_y']
        center_z = config_['center_z']
        size_x = config_['size_x']
        size_y = config_['size_y']
        size_z = config_['size_z']
        covalent_bond_constraints = config_['covalent_bond_constraints']
        
        # Prepare the ligand
        ligand_struct = structure.create_structure_from_smiles(smi)
        ligand_output_file = os.path.join('{task['tmp_run_dir']}', f"ligand_{{ligand_struct.title}}.maegz")
        
        ligprep_settings = LigprepSettings()
        ligprep_settings.set_output_file(ligand_output_file)
        ligprep_job = LigprepJob(ligprep_settings, input_structure=ligand_struct)
        ligprep_job.run()
        ligprep_job.wait()
        
        # Prepare the receptor and ligand structures
        receptor_struct = structure.StructureReader(receptor).next()
        ligand_struct = structure.StructureReader(ligand_output_file).next()
        
        # Set up CovDock settings
        settings = covdock.CovDockSettings()
        settings.set_receptor(receptor_struct)
        settings.set_ligand(ligand_struct)
        
        output_file = '{task['output_path']}'
        settings.set_output_file(output_file)
        settings.set_covalent_bond_atom_pairs(covalent_bond_constraints)
        
        # Specify the ligand binding site as coordinates and box size
        settings.set_site_box_center((center_x, center_y, center_z))
        settings.set_site_box_size((size_x, size_y, size_z))
        
        # Run the CovDock job
        covdock_job = covdock.CovDock(settings)
        covdock_job.run()
        covdock_job.wait()
        
        # Extract the docking scores
        output_structures = list(structure.StructureReader(output_file))
        docking_scores = [struct.property['r_i_docking_score'] for struct in output_structures]
        
        # Save the minimum score to a file
        with open('min_score_covdock.txt', 'w') as score_file:
            score_file.write(str(min(docking_scores)))
        """)
    cmd = ['python3', script_path]
    return cmd

def docking_finish_covdock(task):
    try:
        min_score_file = os.path.join(task['tmp_run_dir'], 'min_score_covdock.txt')
        with open(min_score_file, 'r') as score_file:
            task['score'] = float(score_file.read())

        task['status'] = "success"
    except:
        logging.error("failed parsing")

    return


# Glide SP 
def docking_start_glide_sp(task):
    script_path = os.path.join(task['tmp_run_dir'], 'glide_sp_script.py')
    with open(script_path, 'w') as script_file:
        script_file.write(f"""\
            from schrodinger import structure
            from schrodinger.job import jobcontrol
            from schrodinger.application.glide import glide
            from schrodinger.application.ligprep import LigprepJob, LigprepSettings
            
            config_ = {task['config']}
            smi = '{task['ligand_path']}'
            receptor = config_['receptor']
            center_x = config_['center_x']
            center_y = config_['center_y']
            center_z = config_['center_z']
            size_x = config_['size_x']
            size_y = config_['size_y']
            size_z = config_['size_z']
            
            # Prepare the ligand
            ligand_struct = structure.create_structure_from_smiles(smi)
            ligand_output_file = os.path.join('{task['tmp_run_dir']}', f"ligand_{{ligand_struct.title}}.maegz")
            
            ligprep_settings = LigprepSettings()
            ligprep_settings.set_output_file(ligand_output_file)
            ligprep_job = LigprepJob(ligprep_settings, input_structure=ligand_struct)
            ligprep_job.run()
            ligprep_job.wait()
            
            # Prepare the receptor and ligand structures
            receptor_struct = structure.StructureReader(receptor).next()
            ligand_struct = structure.StructureReader(ligand_output_file).next()
            
            # Set up Glide settings
            settings = glide.GlideSettings()
            settings.set_receptor_file(receptor_struct)
            settings.set_ligand_file(ligand_output_file)
            
            output_file = '{task['output_path']}'
            settings.set_output_file(output_file)
            
            # Specify the ligand binding site as coordinates and box size
            settings.set_site_box_center((center_x, center_y, center_z))
            settings.set_site_box_size((size_x, size_y, size_z))
            
            # Set Glide precision to SP
            settings.set_precision("SP")
            
            # Run the Glide job
            glide_job = glide.Glide(settings)
            glide_job.run()
            glide_job.wait()
            
            # Extract the docking scores
            output_structures = list(structure.StructureReader(output_file))
            docking_scores = [struct.property['r_i_docking_score'] for struct in output_structures]
            
            # Save the minimum score to a file
            with open('min_score.txt', 'w') as score_file:
                score_file.write(str(min(docking_scores)))
            """)
        
        cmd = ['python3', script_path]
        
    return cmd

def docking_finish_glide_sp(task):
    try: 
        min_score_file = os.path.join(task['tmp_run_dir'], 'min_score.txt')
        with open(min_score_file, 'r') as score_file:
            task['score'] = float(score_file.read()) 
        task['status'] = "success"
    except: 
        logging.error("failed parsing")

    return 

# Glide XP 
def docking_start_glide_xp(task):
    script_path = os.path.join(task['tmp_run_dir'], 'glide_xp_script.py')
    with open(script_path, 'w') as script_file:
        script_file.write(f"""\
        from schrodinger import structure
        from schrodinger.job import jobcontrol
        from schrodinger.application.glide import glide
        from schrodinger.application.ligprep import LigprepJob, LigprepSettings
        
        config_ = {task['config']}
        smi = '{task['ligand_path']}'
        receptor = config_['receptor']
        center_x = config_['center_x']
        center_y = config_['center_y']
        center_z = config_['center_z']
        size_x = config_['size_x']
        size_y = config_['size_y']
        size_z = config_['size_z']
        
        # Prepare the ligand
        ligand_struct = structure.create_structure_from_smiles(smi)
        ligand_output_file = os.path.join('{task['tmp_run_dir']}', f"ligand_{{ligand_struct.title}}.maegz")
        
        ligprep_settings = LigprepSettings()
        ligprep_settings.set_output_file(ligand_output_file)
        ligprep_job = LigprepJob(ligprep_settings, input_structure=ligand_struct)
        ligprep_job.run()
        ligprep_job.wait()
        
        # Prepare the receptor and ligand structures
        receptor_struct = structure.StructureReader(receptor).next()
        ligand_struct = structure.StructureReader(ligand_output_file).next()
        
        # Set up Glide settings
        settings = glide.GlideSettings()
        settings.set_receptor_file(receptor_struct)
        settings.set_ligand_file(ligand_output_file)
        
        output_file = '{task['output_path']}'
        settings.set_output_file(output_file)
        
        # Specify the ligand binding site as coordinates and box size
        settings.set_site_box_center((center_x, center_y, center_z))
        settings.set_site_box_size((size_x, size_y, size_z))
        
        # Set Glide precision to XP
        settings.set_precision("XP")
        
        # Run the Glide job
        glide_job = glide.Glide(settings)
        glide_job.run()
        glide_job.wait()
        
        # Extract the docking scores
        output_structures = list(structure.StructureReader(output_file))
        docking_scores = [struct.property['r_i_docking_score'] for struct in output_structures]
        
        # Save the minimum score to a file
        with open('min_score_xp.txt', 'w') as score_file:
            score_file.write(str(min(docking_scores)))
        """)
    cmd = ['python3', script_path]
    return cmd

def docking_finish_glide_xp(task):
    try:
        min_score_file = os.path.join(task['tmp_run_dir'], 'min_score_xp.txt')
        with open(min_score_file, 'r') as score_file:
            task['score'] = float(score_file.read())

        task['status'] = "success"
    except:
        logging.error("failed parsing")

    return


# Glide HTVS
    from schrodinger.job import jobcontrol
    from schrodinger.application.glide import glide
    from schrodinger.application.ligprep import LigprepJob, LigprepSettings

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    smi = task['ligand_path']
    receptor = config_['receptor']
    center_x = config_['center_x']
    center_y = config_['center_y']
    center_z = config_['center_z']
    size_x = config_['size_x']
    size_y = config_['size_y']
    size_z = config_['size_z']

    if receptor.split('.')[-1] != 'maegz':
        logging.error("failed parsing")

    # Prepare the ligand
    ligand_struct = structure.create_structure_from_smiles(smi)
    ligand_output_file = os.path.join(task['tmp_run_dir'], f"ligand_{ligand_struct.title}.maegz")

    ligprep_settings = LigprepSettings()
    ligprep_settings.set_output_file(ligand_output_file)
    ligprep_job = LigprepJob(ligprep_settings, input_structure=ligand_struct)
    ligprep_job.run()
    ligprep_job.wait()

    if ligprep_job.status != jobcontrol.FINISHED:
        logging.error("failed parsing")

    # Prepare the receptor and ligand structures
    receptor_struct = structure.StructureReader(receptor).next()
    ligand_struct = structure.StructureReader(ligand_output_file).next()

    # Set up Glide settings
    settings = glide.GlideSettings()
    settings.set_receptor_file(receptor_struct)
    settings.set_ligand_file(ligand_output_file)

    output_file = task['output_path']
    settings.set_output_file(output_file)

    # Specify the ligand binding site as coordinates and box size
    settings.set_site_box_center((center_x, center_y, center_z))
    settings.set_site_box_size((size_x, size_y, size_z))

    # Set Glide precision to HTVS
    settings.set_precision("HTVS")

    # Run the Glide job
    glide_job = glide.Glide(settings)
    glide_job.run()
    glide_job.wait()

    if glide_job.status != jobcontrol.FINISHED:
        logging.error("failed parsing")

    # Read the output file
    output_structures = list(structure.StructureReader(output_file))

    # Extract the docking scores
    docking_scores = []
    for struct in output_structures:
        docking_score = struct.property['r_i_docking_score']
        docking_scores.append(docking_score)

    task['score'] = min(docking_scores)
    task['status'] = "success"

    return 
def docking_start_glide_htvs(task):
    script_path = os.path.join(task['tmp_run_dir'], 'glide_htvs_script.py')
    with open(script_path, 'w') as script_file:
        script_file.write(f"""\
        from schrodinger import structure
        from schrodinger.job import jobcontrol
        from schrodinger.application.glide import glide
        from schrodinger.application.ligprep import LigprepJob, LigprepSettings
        
        config_ = {task['config']}
        smi = '{task['ligand_path']}'
        receptor = config_['receptor']
        center_x = config_['center_x']
        center_y = config_['center_y']
        center_z = config_['center_z']
        size_x = config_['size_x']
        size_y = config_['size_y']
        size_z = config_['size_z']
        
        # Prepare the ligand
        ligand_struct = structure.create_structure_from_smiles(smi)
        ligand_output_file = os.path.join('{task['tmp_run_dir']}', f"ligand_{{ligand_struct.title}}.maegz")
        
        ligprep_settings = LigprepSettings()
        ligprep_settings.set_output_file(ligand_output_file)
        ligprep_job = LigprepJob(ligprep_settings, input_structure=ligand_struct)
        ligprep_job.run()
        ligprep_job.wait()
        
        # Prepare the receptor and ligand structures
        receptor_struct = structure.StructureReader(receptor).next()
        ligand_struct = structure.StructureReader(ligand_output_file).next()
        
        # Set up Glide settings
        settings = glide.GlideSettings()
        settings.set_receptor_file(receptor_struct)
        settings.set_ligand_file(ligand_output_file)
        
        output_file = '{task['output_path']}'
        settings.set_output_file(output_file)
        
        # Specify the ligand binding site as coordinates and box size
        settings.set_site_box_center((center_x, center_y, center_z))
        settings.set_site_box_size((size_x, size_y, size_z))
        
        # Set Glide precision to HTVS
        settings.set_precision("HTVS")
        
        # Run the Glide job
        glide_job = glide.Glide(settings)
        glide_job.run()
        glide_job.wait()
        
        # Extract the docking scores
        output_structures = list(structure.StructureReader(output_file))
        docking_scores = [struct.property['r_i_docking_score'] for struct in output_structures]
        
        # Save the minimum score to a file
        with open('min_score_htvs.txt', 'w') as score_file:
            score_file.write(str(min(docking_scores)))
        """)
    cmd = ['python3', script_path]
    return cmd

def docking_finish_glide_htvs(task):
    try:
        min_score_file = os.path.join(task['tmp_run_dir'], 'min_score_htvs.txt')
        with open(min_score_file, 'r') as score_file:
            task['score'] = float(score_file.read())

        task['status'] = "success"
    except:
        logging.error("failed parsing")

    return

## PSOvina

def docking_start_PSOVina(task):
    cpus_per_program = str(task['threads_per_docking'])

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = [
            f"{task['tools_path']}/PSOVina",
            '--receptor', config_['receptor'],
            '--ligand', task['ligand_path'],
            '--cpu', cpus_per_program,
            '--exhaustiveness', config_['exhaustiveness'],
            '--center_x', '{}'.format(config_['center_x']),
            '--center_y', '{}'.format(config_['center_y']),
            '--center_z', '{}'.format(config_['center_z']),
            '--size_x',   '{}'.format(config_['size_x']),
            '--size_y',   '{}'.format(config_['size_y']),
            '--size_z',   '{}'.format(config_['size_z']),
            '--out', task['output_path']
         ]
    return cmd

def docking_finish_PSOVina(item, ret):
    match = re.search(r'^\s+1\s+(?P<value>[-0-9.]+)\s+', ret.stdout, flags=re.MULTILINE)
    if(match):
        matches = match.groupdict()
        item['score'] = float(matches['value'])
        item['status'] = "success"
    else:
        item['log']['reason'] = f"Could not find score"
        logging.error(item['log']['reason'])

## *vina

def docking_start_vina(task):
    cpus_per_program = str(task['threads_per_docking'])

    cmd = [
            f"{task['tools_path']}/{task['program']}",
            '--cpu', cpus_per_program,
            '--config', task['config_path'],
            '--ligand', task['ligand_path'],
            '--out', task['output_path']
        ]
    return cmd

def docking_finish_vina(item, ret):
    match = re.search(r'^\s+1\s+(?P<value>[-0-9.]+)\s+', ret.stdout, flags=re.MULTILINE)
    if(match):
        matches = match.groupdict()
        item['score'] = float(matches['value'])
        item['status'] = "success"

        if(item['run_pose_check'] == 1 or item['run_energy_check'] == 1):

            # Load in config file:
            with open(item['config_path']) as fd:
                config_ = dict(read_config_line(line) for line in fd)
            for element in config_:
                if '#' in config_[element]:
                    config_[element] = config_[element].split('#')[0]

            # Convert docking output to sdf (first model only)
            sdf_file = item['output_path'] + '.sdf'
            try:
                ret = subprocess.run(['obabel', '-ipdbqt', item['output_path'], '-f', '1', '-l', '1', '-osdf', '-O', sdf_file],
                                     capture_output=True, text=True, timeout=item['obenergy_timeout'])
            except subprocess.TimeoutExpired as err:
                item['log']['reason'] = 'obabel timed out'
                logging.error(item['log']['reason'])
                logging.error(f"stdout:\n{ret.stdout}\nstderr:{ret.stderr}\n")

            if (item['run_pose_check'] == 1):

                # Prepare receptor file if needed
                pdb_file = os.path.join(item['input_files_dir'], config_['receptor'].replace('pdbqt', 'pdb'))
                if not os.path.exists(pdb_file):
                    logging.error(f"run_pose_check is set to 1, but no {pdb_file} found, exiting...")
                    sys.exit(1)

                # Run pose check
                if(posebusters_check(item, sdf_file, pdb_file, item['pose_check_timeout']) != 1):
                    item['score'] = None
                    item['status'] = "failed"

            if(item['run_energy_check'] == 1):

                # Run energy check
                if(obabel_check_energy(item, sdf_file, item['energy_max']) != 1):
                    item['score'] = None
                    item['status'] = "failed"

    else:
        item['log']['reason'] = f"Could not find score"
        logging.error(item['log']['reason'])

## smina

def docking_start_smina(task):
    cpus_per_program = str(task['threads_per_docking'])
    log_file = os.path.join(task['output_dir'], "out.flexres.pdb")
    atomterms_file = os.path.join(task['output_dir'], "out.atomterms")

    cmd = [
        f"{task['tools_path']}/smina",
        '--cpu', cpus_per_program,
        '--config', task['config_path'],
        '--ligand', task['ligand_path'],
        '--out', task['output_path'],
        '--log', log_file,
        '--atom_terms', atomterms_file
    ]
    return cmd

def docking_finish_smina(item, ret):
    found = 0
    for line in reversed(ret.stdout.splitlines()):
        match = re.search(r'^1\s{4}\s*(?P<value>[-0-9.]+)\s*', line)
        if(match):
            matches = match.groupdict()
            item['score'] = float(matches['value'])
            item['status'] = "success"
            found = 1
            break
    if(found == 0):
        item['log']['reason'] = f"Could not find score"
        logging.error(item['log']['reason'])


## plants

def docking_start_plants(task):

    task['plants_tmp_file'] = os.path.join(task['tmp_run_dir'], "vfvs_tmp.txt")
    shutil.copy(task['config_path'], task['plants_tmp_file'])
        
    with open(task['plants_tmp_file'], 'a+') as f:
        f.writelines('ligand_file {}\n'.format(task['ligand_path']))
        f.writelines('output_dir {}\n'.format(task['output_path']))


    cmd = ['{}/PLANTS'.format(task['tools_path']),
            '--mode', 'screen',
            task['plants_tmp_file']
    ]

    return cmd

def docking_finish_plants(item, ret):

    try:
        plants_cmd = ret.stdout.split('\n')
        plants_cmd = [x for x in plants_cmd if 'best score:' in x][-1]
        item['score'] = float(plants_cmd.split(' ')[-1])
        item['status'] = "success"
    except:
        logging.error("failed parsing")

## adfr

def docking_start_adfr(item):
    with open(item['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    cmd = ['{}/adfr'.format(item['tools_path']),
           '-t', '{}'.format(config_['receptor']),
           '-l', '{}'.format(item['ligand_path']),
           '--jobName', '{}'.format(item['output_path'])
           ]
    return cmd

def docking_finish_adfr(item, ret):

    try:
        docking_out = ret.stdout
        docking_scores = []
        for line_item in docking_out:
            A = line_item.split(' ')
            A = [x for x in A if x != '']
            try:
                _, a_2, _ = float(A[0]), float(A[1]), float(A[2])
            except:
                continue
            docking_scores.append(float(a_2))

        item['score'] = min(docking_scores)
        item['status'] = "success"
    except:
        logging.error("failed parsing")

## AutodockZN

def docking_start_autodockzn(task):

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/AutodockVina_1.2'.format(task['tools_path']),
           '--ligand', '{}'.format(task['ligand_path']),
           '--maps', config_['afinit_maps_name'],
           '--scoring', 'ad4',
           '--exhaustiveness', '{}'.format(config_['exhaustiveness']),
           '--out', '{}'.format(task['output_path'])]
    return cmd

## gnina

def docking_start_gnina(task):

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/gnina'.format(task['tools_path']),
               '-r', config_['receptor'],
               '-l', '{}'.format(task['ligand_path']),
               '--exhaustiveness', '{}'.format(config_['exhaustiveness']),
               '--center_x', '{}'.format(config_['center_x']),
               '--center_y', '{}'.format(config_['center_y']),
               '--center_z', '{}'.format(config_['center_z']),
               '--size_x',   '{}'.format(config_['size_x']),
               '--size_y',   '{}'.format(config_['size_y']),
               '--size_z',   '{}'.format(config_['size_z']),
               '--out', '{}'.format(task['output_path'])]
    return cmd

## rDock

def docking_start_rdock(task):
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = [ 'rbdock',
                '-i', task['ligand_path'],
                '-o', task['output_path'],
                '-r', config_['rdock_config'],
                '-p', config_['dock_prm'],
                '-n', config_['runs']]
    return cmd

def docking_finish_rdock(item, ret):

    try:
        with open(item['output_path'], 'r') as f:
            lines = f.readlines()
        score = []
        for i, item in enumerate(lines):
            if item.strip() == '>  <SCORE>':
                score.append(float(lines[i+1]))
        item['score'] = min(score)
        item['status'] = "success"
    except:
        logging.error("failed parsing")


## M-Dock

def docking_start_mdock(task):

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/MDock_Linux'.format(task['tools_path']),
            config_['protein_name'],
            task['ligand_path'],
           '-param', config_['mdock_config']
          ]

    return cmd


def docking_finish_mdock(item, ret):

    try:
        docking_scores = []

        output_file = os.path.join(item['tmp_run_dir_input'], "mdock_dock.out")
        with open(output_file, 'r') as f:
            lines = f.readlines()

        for item in lines:
            docking_scores.append( float([x for x in item.split(' ') if x != ''][4]))

        shutil.move(output_file, item['output_dir'])
        mol_output_file = os.path.join(item['tmp_run_dir_input'], "mdock_dock.mol2")
        shutil.move(mol_output_file, item['output_path'])

        item['score'] = min(docking_scores)
        item['status'] = "success"
    except:
        logging.error("failed parsing")

## MCDock

def docking_start_mcdock(task):

    with open(task['config_path']) as fd:
            config_ = dict(read_config_line(line) for line in fd)

    cmd = ['{}/mcdock'.format(task['tools_path']),
        '--target', config_['protein_name'],
        '--ligand', task['ligand_path']]

    return cmd

def docking_finish_mcdock(item, ret):

    try:
        output_file = os.path.join(item['tmp_run_dir_input'], "out.xyz")

        with open(output_file, 'r') as f:
            lines = f.readlines()

        lines = [x for x in lines if 'Binding Energy' in x]
        binding_energies = []
        for item in lines:
            binding_energies.append(float(item.split(' ')[2].split('\t')[0]))

        item['score'] = min(binding_energies)
        item['status'] = "success"

        shutil.move(output_file, item['output_path'])
    except:
        logging.error("failed parsing")

## LigandFit

def docking_start_ligandfit(task):

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/ligandfit'.format(task['tools_path']),
           'data=', config_['receptor_mtz'],
           'model=', config_['receptor'],
           'ligand', task['ligand_path'],
           'search_center=', config_['center_x'], config_['center_y'], config_['center_z']]

    return cmd

def docking_finish_ligandfit(item, ret):

    run_pdb = os.path.join(item['tmp_run_dir_input'], "LigandFit_run_1_", "ligand_fit_1.pdb")
    run_log = os.path.join(item['tmp_run_dir_input'], "LigandFit_run_1_", "ligand_1_1.log")

    try:
        with open(run_log, 'r') as f:
            lines = f.readlines()
        lines = [x for x in lines if 'Best score' in x]
        scores = []
        for item in lines:
            scores.append( float([x for x in item.split(' ') if x != ''][-2]) )

        item['score'] = min(scores)
        item['status'] = "success"

        shutil.move(run_pdb, item['output_path'])
        shutil.move(run_log, item['output_dir'])
    except:
        logging.error("failed parsing")

    # TODO
    os.system('rm -rf LigandFit_run_1_')

## ledock

def docking_start_ledock(item):

    item['ledock_tmp_file'] = os.path.join(item['tmp_run_dir'], "vfvs_tmp.in")
    item['ledock_tmp_file_list'] = os.path.join(item['tmp_run_dir'], "vfvs_tmp.list")

    cmd = [
        '{}/ledock'.format(item['tools_path']),
        '{}'.format(item['ledock_tmp_file'])
    ]

    with open(item['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]: config_[item] = config_[item].split('#')[0]

    docking_file = cmd[-1]
    # ligand_list_file = docking_file.split('.')[0] + '.list'

    with open(docking_file, 'w') as f:
        f.writelines(['Receptor'])
        f.writelines([config_['receptor'] + '\n'])
        f.writelines(['RMSD'])
        f.writelines([config_['rmsd'] + '\n'])
        f.writelines(['Binding pocket'])
        f.writelines(['{} {}'.format( config_['min_x'], config_['max_x']) ])
        f.writelines(['{} {}'.format( config_['min_y'], config_['max_y']) ])
        f.writelines(['{} {}\n'.format( config_['min_z'], config_['max_z']) ])
        f.writelines(['Number of binding poses'])
        f.writelines([config_['n_poses'] + '\n'])
        f.writelines(['Ligands list'])
        f.writelines([item['ledock_tmp_file_list'] + '\n'])
        f.writelines(['END'])

    with open(item['ledock_tmp_file_list'], 'w') as f:
        f.writelines(item['ligand_path'])

    return cmd

def docking_finish_ledock(item, ret):

    try:
        ligand_filename = item['ligand_path'].split('/')[-1]
        ligand_base = ligand_filename.split('.')[0]

        run_dok = os.path.join(item['tmp_run_dir_input'], "ligands", f"{ligand_base}.dok")

        with open(run_dok, 'r') as f:
            lines = f.readlines()
        lines = [x for x in lines if 'Score' in x]
        scores = []
        for item in lines:
            A = item.split('Score')[-1].strip().split(': ')[1].split(' ')[0]
            scores.append(float(A))
        item['score'] = min(scores)
        item['status'] = "success"

        shutil.move(run_dok, item['output_dir'])

    except:
        logging.error("failed parsing")


## gold

def docking_start_gold(item):

    item['gold_tmp_file'] = os.path.join(item['tmp_run_dir'], "vfvs_tmp.conf")
    item['gold_tmp_dir'] = os.path.join(item['tmp_run_dir'], "vfvs_tmp")

    with open(item['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]: config_[item] = config_[item].split('#')[0]

    with open(item['gold_tmp_file'], 'w') as f:
        f.writelines(['  GOLD CONFIGURATION FILE\n'])
        f.writelines(['  AUTOMATIC SETTINGS'])
        f.writelines(['autoscale = 1\n'])
        f.writelines(['  POPULATION'])
        f.writelines(['popsiz = auto'])
        f.writelines(['select_pressure = auto'])
        f.writelines(['n_islands = auto'])
        f.writelines(['maxops = auto'])
        f.writelines(['niche_siz = auto\n'])
        f.writelines(['  GENETIC OPERATORS'])
        f.writelines(['pt_crosswt = auto'])
        f.writelines(['allele_mutatewt = auto'])
        f.writelines(['migratewt = auto\n'])
        f.writelines(['  FLOOD FILL'])
        f.writelines(['radius = {}'.format(config_['radius'])])
        f.writelines(['origin = {}   {}   {}'.format(config_['center_x'], config_['center_y'], config_['center_z'])])
        f.writelines(['do_cavity = 0'])
        f.writelines(['floodfill_center = point\n'])
        f.writelines(['   DATA FILES'])
        f.writelines(['ligand_data_file {} 10'.format(item['ligand_path'])])
        f.writelines(['param_file = DEFAULT'])
        f.writelines(['set_ligand_atom_types = 1'])
        f.writelines(['set_protein_atom_types = 0'])
        f.writelines(['directory = {}'.format(item['gold_tmp_dir'])])
        f.writelines(['tordist_file = DEFAULT'])
        f.writelines(['make_subdirs = 0'])
        f.writelines(['save_lone_pairs = 1'])
        f.writelines(['fit_points_file = fit_pts.mol2'])
        f.writelines(['read_fitpts = 0'])
        f.writelines(['bestranking_list_filename = bestranking.lst\n'])
        f.writelines(['   FLAGS'])
        f.writelines(['internal_ligand_h_bonds = 1'])
        f.writelines(['flip_free_corners = 1'])
        f.writelines(['match_ring_templates = 1'])
        f.writelines(['flip_amide_bonds = 0'])
        f.writelines(['flip_planar_n = 1 flip_ring_NRR flip_ring_NHR'])
        f.writelines(['flip_pyramidal_n = 0'])
        f.writelines(['rotate_carboxylic_oh = flip'])
        f.writelines(['use_tordist = 1'])
        f.writelines(['postprocess_bonds = 1'])
        f.writelines(['rotatable_bond_override_file = DEFAULT'])
        f.writelines(['solvate_all = 1\n'])
        f.writelines(['   TERMINATION'])
        f.writelines(['early_termination = 1'])
        f.writelines(['n_top_solutions = 3'])
        f.writelines(['rms_tolerance = 1.5\n'])
        f.writelines(['   CONSTRAINTS'])
        f.writelines(['force_constraints = 0\n'])
        f.writelines(['   COVALENT BONDING'])
        f.writelines(['covalent = 0\n'])
        f.writelines(['   SAVE OPTIONS'])
        f.writelines(['save_score_in_file = 1'])
        f.writelines(['save_protein_torsions = 1\n'])
        f.writelines(['  FITNESS FUNCTION SETTINGS'])
        f.writelines(['initial_virtual_pt_match_max = 4'])
        f.writelines(['relative_ligand_energy = 1'])
        f.writelines(['gold_fitfunc_path = goldscore'])
        f.writelines(['score_param_file = DEFAULT\n'])
        f.writelines(['  PROTEIN DATA'])
        f.writelines(['protein_datafile = {}'.format(config_['receptor'])])


    cmd = ['{}/gold_auto'.format(item['tools_path']), '{}'.format(item['gold_tmp_file'])]

    return cmd


def docking_finish_gold(item, ret):
    try:
        # TODO -- fix all of this...

        run_output = os.path.join(item['gold_tmp_dir'], "ligand_m1.rnk")
        run_pose = os.path.join(item['gold_tmp_dir'], "gold_ligand_m1.mol2")

        with open(run_output, 'r') as f:
            lines = f.readlines()
            docking_score = float([x for x in lines[-1].split(' ') if x!=''][1])
            item['score'] = min(docking_score)
            item['status'] = "success"

        shutil.move(run_pose, item['output_dir'])
        shutil.move(run_output, item['output_dir'])

    except:
        logging.error("failed parsing")


## iGemDock

def docking_start_igemdock(task):

    task['igemdock_temp_dir'] = os.path.join(task['tmp_run_dir'], "vfvs_tmp")

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/mod_ga'.format(task['tools_path']),
           config_['exhaustiveness'],
           config_['receptor'],
           task['ligand_path'],
           '-d', item['igemdock_temp_dir']
    ]

    return cmd


def docking_finish_igemdock(task, ret):

    try:

        docked_pose = os.listdir(os.path.join(task['igemdock_temp_dir'], ''))[0]
        with open(docked_pose, 'r') as f:
            lines = f.readlines()

        docking_score = lines[4]
        docking_score = float([x for x in docking_score.split(' ') if x != ''][1])

        shutil.move(docked_pose, task['output_path'])

        task['score'] = min(docking_score)
        task['status'] = "success"
    except:
        logging.error("failed parsing")

## idock

def docking_start_idock(item):

    with open(item['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['{}/idock'.format(item['tools_path']),
           '--receptor', config_['receptor'],
           '--ligand', item['ligand_path'],
           '--center_x', config_['center_x'],
           '--center_y', config_['center_y'],
           '--center_z', config_['center_z'],
           '--size_x', config_['size_x'],
           '--size_y', config_['size_y'],
           '--size_z', config_['size_z'],
           '--out', '{}'.format(item['output_path'])]

    return cmd

def docking_finish_idock(item, ret):
    try:
        docking_out = ret.stdout
        docking_out = float([x for x in docking_out.split(' ') if x != ''][-2])
        item['score'] = min(docking_out)
        item['status'] = "success"
    except:
        logging.error("failed parsing")

## GalaxyDock3

def docking_start_galaxydock3(task):

    task['galaxydock3_tmp_file'] = os.path.join(task['tmp_run_dir'], "vfvs_tmp.in")
    task['ligdock_prefix'] = "vfvs_tmp"

    cmd = [
            '{}/GalaxyDock3'.format(task['tools_path']),
            task['galaxydock3_tmp_file']
    ]

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]: config_[item] = config_[item].split('#')[0]

    with open(cmd[-1], 'w') as f:
        f.writelines(['!=============================================='])
        f.writelines(['! I/O Parameters'])
        f.writelines(['!=============================================='])
        f.writelines(['data_directory    ./'])
        f.writelines(['infile_pdb        {}'.format(config_['receptor'])])
        f.writelines(['infile_ligand        {}'.format(task['ligand_path'])])
        f.writelines(['top_type          polarh'])
        f.writelines(['fix_type          all'])
        f.writelines(['ligdock_prefix    {}'.format(task['ligdock_prefix'])])
        f.writelines(['!=============================================='])
        f.writelines(['! Grid Options'])
        f.writelines(['!=============================================='])
        f.writelines(['grid_box_cntr     {} {} {}'.format(config_['grid_box_cntr'].split(' ')[0], config_['grid_box_cntr'].split(' ')[1], config_['grid_box_cntr'].split(' ')[2])])
        f.writelines(['grid_n_elem       {} {} {}'.format(config_['grid_n_elem'].split(' ')[0], config_['grid_n_elem'].split(' ')[1], config_['grid_n_elem'].split(' ')[2])])
        f.writelines(['grid_width        {}'.format(config_['grid_width'])])
        f.writelines(['!=============================================='])
        f.writelines(['! Energy Parameters'])
        f.writelines(['!=============================================='])
        f.writelines(['weight_type              GalaxyDock3'])
        f.writelines(['!=============================================='])
        f.writelines(['! Initial Bank Parameters'])
        f.writelines(['!=============================================='])
        f.writelines(['first_bank               rand'])
        f.writelines(['max_trial                {}'.format(config_['max_trial'])])
        f.writelines(['e0max                    1000.0'])
        f.writelines(['e1max                    1000000.0'])
        f.writelines(['n_proc 1'])


    return cmd

def docking_finish_galaxydock3(item, ret):
    try:

        info_file = os.path.join(item['tmp_run_dir_input'], f"{item['ligdock_prefix']}_fb.E.info")
        mol_file = os.path.join(item['tmp_run_dir_input'], f"{item['ligdock_prefix']}_fb.mol2")

        with open(info_file, 'r') as f:
            lines = f.readlines()
        lines = lines[3: ]
        docking_scores = []
        for item in lines:
            try:
                A = item.split(' ')
                A = [x for x in A if x != '']
                docking_scores.append(float(A[5]))
            except:
                continue

        shutil.move(info_file, item['output_dir'])
        shutil.move(mol_file, item['output_dir'])

        item['score'] = min(docking_scores)
        item['status'] = "success"
    except:
        logging.error("failed parsing")


## Autodock

def docking_start_autodock_gpu(task):
    return docking_start_autodock(task, "gpu")

def docking_start_autodock_cpu(task):
    return docking_start_autodock(task, "cpu")

def docking_start_autodock(item, arch_type):

    with open(item['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]: config_[item] = config_[item].split('#')[0]

    cmd = ['{}/autodock_{}'.format(item['tools_path'], arch_type),
           '--ffile', config_['receptor'],
           '--lfile', item['ligand_path']]

    return cmd

def docking_finish_autodock(item, ret):
    try :
        output = ret.stdout.split('\n')[-6]
        lines = [x.strip() for x in output if 'best energy' in x][0]
        docking_score = float(lines.split(',')[1].split(' ')[-2])
        item['score'] = min(docking_score)
        item['status'] = "success"
    except:
        logging.error("failed parsing")

## FRED

def docking_start_fred(task):
    if not os.path.exists('./oe_license.txt'): 
        raise Exception('OpenEye license file (oe_license.txt) not found in tools path')

    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)

    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    cmd = ['python {}/dock_fred.py'.format(task['tools_path']),
           '--receptor-fn', config_['receptor'],
           '--ligand-fn', task['ligand_path'],
           '--center_x', config_['center_x'],
           '--center_y', config_['center_y'],
           '--center_z', config_['center_z'],
           '--radius', max([config_['size_x'], config_['size_y'], config_['size_z']]),
           '--num-poses', config_['exhaustiveness'],
           '--out', '{}'.format(task['output_path'])]
        
    return cmd

def docking_finish_fred(task, ret):
    cmd_scoring = scoring_start_vina(task)
    try: 
        subprocess.run(cmd_scoring, capture_output=True,check=True)
        scoring_finish_vina(item=task, ret=ret)
    except:
        logging.error("failed scoring")

    return


# Scoring functions: 
def convert_ligand_format(ligand_, new_format): 
    """Converts a ligand file to a different file format using the Open Babel tool.

        Args:
            ligand_ (str): The path to the input ligand file.
            new_format (str): The desired output format for the ligand file.
    
        Returns:
            None
    
        Raises:
            Exception: If the input file does not exist, or if the Open Babel tool is not installed.
    
        Examples:
            To convert a ligand file from mol2 format to pdbqt format:
            >>> convert_ligand_format('./ligands/ligand1.mol2', 'pdbqt')
    """
    input_format = ligand_.split('.')[-1]
    os.system('obabel {} -O {}'.format(ligand_, ligand_.replace(input_format, new_format)))

    
## nnscore2.0
def scoring_start_nnscore2(task): 
    
    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format =  task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format( task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    vina_loc = '{}/vina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('export VINA_EXEC={}; python {}/NNScore2.py -receptor {} -ligand {} -vina_executable $VINA_EXEC > output.txt'.format(vina_loc, task['tools_path'], config_['receptor'], task['output_path']))
    
    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_nnscore2(item, ret): 
    
    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        scores = [x for x in lines if 'Best Score:' in x]
        scores = [A.split('(')[-1].split(')')[0] for A in scores]
        item['score'] = min(scores)   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
## GOLD scoring
def _scoring_start_gold(task, scoring_function: str):
    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format =  task['output_path'].split('.')[-1]
    if lig_format not in ['mol2', 'mol', 'mdl', 'sdf']: 
        convert_ligand_format(task['output_path'], 'mol2')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'mol2') 

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    input_conf = os.path.join(task['tmp_run_dir'], 'input.conf')
    gold_loc = '{}/gold_auto'.format(task['tools_path'])

    output_dir = 'gold_output'
    os.mkdir(os.path.join(task['tmp_run_dir'], output_dir))
    with open(input_conf, mode='w') as f:
        f.writelines([
            'protein_datafile = {}\n'.format(config_['receptor']),
            'ligand_data_file = {} 10\n'.format(task['output_path']),
            'param_file = DEFAULT\n',
            f'directory = {output_dir}\n',
            f'gold_fitfunc_path {scoring_function}\n',
            'run_flag = RESCORE\n',
        ])
    
    with open(run_sh_script, 'w') as f:
        f.writelines(f'{gold_loc} {input_conf}')

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 

    return cmd 

def _scoring_finish_gold(item, ret):
    try:
        with open(f"{item['tmp_run_dir']}/gold_output/rescore.log") as rescore_log:
            for line in rescore_log:
                pass
        last_line = [i for i in line.split(sep=' ') if i]
        item['score'] = last_line[4]
        item['status'] = "success"
    except:
        logging.error("failed parsing")

def scoring_start_gold_asp(task):
    return _scoring_start_gold(task=task, scoring_function='asp')

def scoring_finish_gold_asp(item, ret):
    return _scoring_finish_gold(item=item, ret=ret)

def scoring_start_gold_chemscore(task):
    return _scoring_start_gold(task=task, scoring_function='chemscore')

def scoring_finish_gold_chemscore(item, ret):
    return _scoring_finish_gold(item=item, ret=ret)

def scoring_start_gold_goldscore(task):
    return _scoring_start_gold(task=task, scoring_function='goldscore')

def scoring_finish_gold_goldscore(item, ret):
    return _scoring_finish_gold(item=item, ret=ret)

def scoring_start_gold_plp(task):
    return _scoring_start_gold(task=task, scoring_function='plp')

def scoring_finish_gold_plp(item, ret):
    return _scoring_finish_gold(item=item, ret=ret)

## rf-score-vs

def scoring_start_rf(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    rf_score_vs_loc = '{}/rf-score-vs'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} {} -O {}/ligands_rescored.pdbqt'.format(rf_score_vs_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))
        f.writelines('{} --receptor {} {} -ocsv > {}/temp.csv'.format(rf_score_vs_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_rf(item, ret): 

    try:    
        with open('{}/temp.csv'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        rf_scores = []
        for line_item in lines[1: ]: 
            rf_scores.append( float(line_item.split(',')[-1]) )
        item['score'] = min(rf_scores)   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

## smina scoring

def scoring_start_smina(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    smina_loc = '{}/smina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} -l {} --score_only > {}/output.txt'.format(smina_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_smina(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        smina_score = float([x for x in lines if 'Affinity' in x][0].split(' ')[1])
        item['score'] = smina_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

## gnina scoring

def scoring_start_gnina(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    smina_loc = '{}/gnina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} -l {} --score_only > {}/output.txt'.format(smina_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_gnina(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        smina_score = float([x for x in lines if 'Affinity' in x][0].split(' ')[1])
        item['score'] = smina_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

## ad4 scoring
def scoring_start_ad4(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    smina_loc = '{}/smina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} -l {} --score_only --scoring ad4_scoring > {}/output.txt'.format(smina_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_ad4(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        ad4_score = float([x for x in lines if 'Affinity' in x][0].split(' ')[1])
        item['score'] = ad4_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")


# vinardo scoring
def scoring_start_vinardo(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
            
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    smina_loc = '{}/smina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} -l {} --score_only --scoring vinardo > {}/output.txt'.format(smina_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_vinardo(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        vinardo_score = float([x for x in lines if 'Affinity' in x][0].split(' ')[1])
        item['score'] = vinardo_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

# vina scoring 
def scoring_start_vina(task):

    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]
    
    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'pdbqt': 
        convert_ligand_format(task['output_path'], 'pdbqt')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'pdbqt')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    smina_loc = '{}/smina'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --receptor {} -l {} --score_only --scoring vina > {}/output.txt'.format(smina_loc, config_['receptor'], task['output_path'], task['tmp_run_dir']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_vina(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        vina_score = float([x for x in lines if 'Affinity' in x][0].split(' ')[1])
        item['score'] = vina_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

# PLANTS chemplp
def scoring_start_PLANTS_chemplp(task):
    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'mol2': 
        convert_ligand_format(task['output_path'], 'mol2')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    plants_loc = '{}/PLANTS'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --mode rescore --config_file {}/plants_config > {}/output.txt'.format(plants_loc, task['tmp_run_dir'], task['tmp_run_dir']))
        
    with open(os.path.join(task['tmp_run_dir'], "plants_config"), 'w') as f:
        f.writelines('scoring_function         chemplp\n')
        f.writelines('protein_file             {}\n'.format(config_['receptor']))
        f.writelines('ligand_file              {}\n'.format(task['output_path']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_PLANTS_chemplp(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        chemplp_score = float([x for x in lines if 'best score:' in x][-1].split(' ')[-1])
        item['score'] = chemplp_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

# PLANTS plp
def scoring_start_PLANTS_plp(task):
    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'mol2': 
        convert_ligand_format(task['output_path'], 'mol2')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    plants_loc = '{}/PLANTS'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --mode rescore --config_file {}/plants_config > {}/output.txt'.format(plants_loc, task['tmp_run_dir'], task['tmp_run_dir']))
        
    with open(os.path.join(task['tmp_run_dir'], "plants_config"), 'w') as f:
        f.writelines('scoring_function         plp\n')
        f.writelines('protein_file             {}\n'.format(config_['receptor']))
        f.writelines('ligand_file              {}\n'.format(task['output_path']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_PLANTS_plp(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        plp_score = float([x for x in lines if 'best score:' in x][-1].split(' ')[-1])
        item['score'] = plp_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")
        
# PLANTS plp95
def scoring_start_PLANTS_plp95(task):
    # Load in config file: 
    with open(task['config_path']) as fd:
        config_ = dict(read_config_line(line) for line in fd)
    for item in config_:
        if '#' in config_[item]:
            config_[item] = config_[item].split('#')[0]

    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'mol2': 
        convert_ligand_format(task['output_path'], 'mol2')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    plants_loc = '{}/PLANTS'.format(task['tools_path'])

    with open(run_sh_script, 'w') as f:        
        f.writelines('{} --mode rescore --config_file {}/plants_config > {}/output.txt'.format(plants_loc, task['tmp_run_dir'], task['tmp_run_dir']))
        
    with open(os.path.join(task['tmp_run_dir'], "plants_config"), 'w') as f:
        f.writelines('scoring_function         plp95\n')
        f.writelines('protein_file             {}\n'.format(config_['receptor']))
        f.writelines('ligand_file              {}\n'.format(task['output_path']))

    os.system('chmod 0700 {}'.format(run_sh_script))
    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_PLANTS_plp95(item, ret): 

    try:    
        with open('{}/output.txt'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        plp95_score = float([x for x in lines if 'best score:' in x][-1].split(' ')[-1])
        item['score'] = plp95_score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")


# Dock6 Contact Score
def scoring_start_dock6_contact_score(task):
    # Load in config file: 
    with open(task['config_path']) as fd:
        config = dict(read_config_line(line) for line in fd)
    for item in config:
        if '#' in config[item]:
            config[item] = config[item].split('#')[0]

    # Convert ligand format if needed:
    lig_format = task['output_path'].split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(task['output_path'], 'mol2')
        task['output_path'] = task['output_path'].replace(task['output_path'], 'mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")

    with open(run_sh_script, 'w') as f:        
        f.writelines(['export Chimera={}\n'.format(config['chimera_path'])])
        f.writelines(['export DOCK6={}\n'.format(config['dock6_path'])])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(config['receptor'])])

        # Generate INSPH, box.in, grid.in, and Contact_Score.in files
        insph_path = os.path.join(task['tmp_run_dir'], 'INSPH')
        with open(insph_path, 'w') as f:
            f.writelines('rec.ms\n')
            f.writelines('R\n')
            f.writelines('X\n')
            f.writelines('0.0\n')
            f.writelines('4.0\n')
            f.writelines('1.4\n')
            f.writelines('rec.sph\n')

        box_path = os.path.join(task['tmp_run_dir'], 'box.in')
        with open(box_path, 'w') as f:
            f.writelines('N\n')
            f.writelines('U\n')
            f.writelines('{}   {}    {}\n'.format(config['center_x'], config['center_y'], config['center_z']))
            f.writelines('{} {} {}\n'.format(config['size_x'], config['size_y'], config['size_z']))
            f.writelines('rec_box.pdb\n')

        grid_path = os.path.join(task['tmp_run_dir'], 'grid.in')
        with open(grid_path, 'w') as f:
            f.writelines('compute_grids                  yes\n')
            f.writelines('energy_score                   yes\n')
            f.writelines('energy_cutoff_distance         9999\n')
            f.writelines('atom_model                     a\n')
            f.writelines('bump_filter                    yes\n')
            f.writelines('receptor_file                  {}\n'.format(config['receptor']))
            f.writelines('box_file                       rec_box.pdb\n')
            f.writelines('vdw_definition_file            {}/parameters/vdw_AMBER_parm99.defn\n'.format(config['dock6_path']))
            f.writelines('score_grid_prefix              grid\n')
            f.writelines('grid_spacing                   0.3\n')
            f.writelines('output_molecule                no\n')
            f.writelines('contact_score                  yes\n')
            f.writelines('attractive_exponent            6\n')
            f.writelines('repulsive_exponent             12\n')
            f.writelines('distance_dielectric            yes\n')
            f.writelines('dielectric_factor              4\n')
            f.writelines('bump_overlap                   0.75\n')
            f.writelines('contact_cutoff_distance        4.5\n')
            
        contact_score_path = os.path.join(task['tmp_run_dir'], 'Contact_Score.in')
        with open(contact_score_path, 'w') as f:
            f.writelines(['conformer_search_type                                        rigid\n'])
            f.writelines(['use_internal_energy                                          yes\n'])
            f.writelines(['internal_energy_rep_exp                                      12\n'])
            f.writelines(['internal_energy_cutoff                                       100.0\n'])
            f.writelines(['ligand_atom_file                                             {}\n'.format(task['output_path'])])
            f.writelines(['limit_max_ligands                                            no\n'])
            f.writelines(['skip_molecule                                                no\n'])
            f.writelines(['read_mol_solvation                                           no\n'])
            f.writelines(['calculate_rmsd                                               no\n'])
            f.writelines(['use_database_filter                                          no\n'])
            f.writelines(['orient_ligand                                                no\n'])
            f.writelines(['bump_filter                                                  no\n'])
            f.writelines(['score_molecules                                              yes\n'])
            f.writelines(['contact_score_primary                                        yes\n'])
            f.writelines(['contact_score_secondary                                      no\n'])
            f.writelines(['contact_score_cutoff_distance                                4.5\n'])
            f.writelines(['contact_score_clash_overlap                                  0.75\n'])
            f.writelines(['contact_score_clash_penalty                                  50\n'])
            f.writelines(['contact_score_grid_prefix                                    grid\n'])
            f.writelines(['grid_score_secondary                                         no\n'])
            f.writelines(['multigrid_score_secondary                                    no\n'])
            f.writelines(['dock3.5_score_secondary                                      no\n'])
            f.writelines(['continuous_score_secondary                                   no\n'])
            f.writelines(['footprint_similarity_score_secondary                         no\n'])
            f.writelines(['pharmacophore_score_secondary                                no\n'])
            f.writelines(['descriptor_score_secondary                                   no\n'])
            f.writelines(['gbsa_zou_score_secondary                                     no\n'])
            f.writelines(['gbsa_hawkins_score_secondary                                 no\n'])
            f.writelines(['SASA_score_secondary                                         no\n'])
            f.writelines(['amber_score_secondary                                        no\n'])
            f.writelines(['minimize_ligand                                              yes\n'])
            f.writelines(['simplex_max_iterations                                       1000\n'])
            f.writelines(['simplex_tors_premin_iterations                               0\n'])
            f.writelines(['simplex_max_cycles                                           1\n'])
            f.writelines(['simplex_score_converge                                       0.1\n'])
            f.writelines(['simplex_cycle_converge                                       1.0\n'])
            f.writelines(['simplex_trans_step                                           1.0\n'])
            f.writelines(['simplex_rot_step                                             0.1\n'])
            f.writelines(['simplex_tors_step                                            10.0\n'])
            f.writelines(['simplex_random_seed                                          0\n'])
            f.writelines(['simplex_restraint_min                                        no\n'])
            f.writelines(['atom_model                                                   all\n'])
            f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(config['dock6_path'])])
            f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(config['dock6_path'])])
            f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(config['dock6_path'])])
            f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
            f.writelines(['write_orientations                                           no\n'])
            f.writelines(['num_scored_conformers                                        1\n'])
            f.writelines(['rank_ligands                                                 no\n'])

    os.system('chmod 0700 {}'.format(run_sh_script))

    cmd = ['./{}'.format(run_sh_script)] 
    
    return cmd 

def scoring_finish_dock6_contact_score(item, ret): 
    try:    
        with open('{}/ligand_out_scored.mol2'.format(item['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()
        score = float([x for x in lines[2].split(' ') if x != ''][-1])
        item['score'] = score   
        item['status'] = "success"
    except: 
        logging.error("failed parsing")

def scoring_start_dock6_continuous_score(task):
    config = load_config(config_path=task['config_path'])
    task['output_path'] = format_ligand(ligand_path=task['output_path'], file_format='mol2')
    
    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    with open(run_sh_script, 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(config['chimera_path'])])
        f.writelines(['export DOCK6={}\n'.format(config['dock6_path'])])
        f.writelines(['$Chimera/bin/chimera --nogui {} {}/dockprep.py\n'.format(config['receptor'], task['tools_path'])])
        f.writelines(['$DOCK6/bin/dock6 -i Continuous_Score.in\n'])

    continuous_score_path = os.path.join(task['tmp_run_dir'], 'Continuous_Score.in')    
    with open(continuous_score_path, 'w') as f: 
        f.writelines(['conformer_search_type                                        rigid\n'])
        f.writelines(['use_internal_energy                                          yes\n'])
        f.writelines(['internal_energy_rep_exp                                      12\n'])
        f.writelines(['internal_energy_cutoff                                       100.0\n'])
        f.writelines(['ligand_atom_file                                             {}\n'.format(task['output_path'])])
        f.writelines(['limit_max_ligands                                            no\n'])
        f.writelines(['skip_molecule                                                no\n'])
        f.writelines(['read_mol_solvation                                           no\n'])
        f.writelines(['calculate_rmsd                                               no\n'])
        f.writelines(['use_database_filter                                          no\n'])
        f.writelines(['orient_ligand                                                no\n'])
        f.writelines(['bump_filter                                                  no\n'])
        f.writelines(['score_molecules                                              yes\n'])
        f.writelines(['contact_score_primary                                        no\n'])
        f.writelines(['contact_score_secondary                                      no\n'])
        f.writelines(['grid_score_primary                                           no\n'])
        f.writelines(['grid_score_secondary                                         no\n'])
        f.writelines(['multigrid_score_primary                                      no\n'])
        f.writelines(['multigrid_score_secondary                                    no\n'])
        f.writelines(['dock3.5_score_primary                                        no\n'])
        f.writelines(['dock3.5_score_secondary                                      no\n'])
        f.writelines(['continuous_score_primary                                     yes\n'])
        f.writelines(['continuous_score_secondary                                   no\n'])
        f.writelines(['cont_score_rec_filename                                      rec_charged.mol2\n'])
        f.writelines(['cont_score_att_exp                                           6\n'])
        f.writelines(['cont_score_rep_exp                                           12\n'])
        f.writelines(['cont_score_rep_rad_scale                                     1.0\n'])
        f.writelines(['cont_score_use_dist_dep_dielectric                           yes\n'])
        f.writelines(['cont_score_dielectric                                        4.0\n'])
        f.writelines(['cont_score_vdw_scale                                         yes\n'])
        f.writelines(['cont_score_turn_off_vdw                                      yes\n'])
        f.writelines(['cont_score_es_scale                                          1.0\n'])
        f.writelines(['footprint_similarity_score_secondary                         no\n'])
        f.writelines(['pharmacophore_score_secondary                                no\n'])
        f.writelines(['descriptor_score_secondary                                   no\n'])
        f.writelines(['gbsa_zou_score_secondary                                     no\n'])
        f.writelines(['gbsa_hawkins_score_secondary                                 no\n'])
        f.writelines(['SASA_score_secondary                                         no\n'])
        f.writelines(['amber_score_secondary                                        no\n'])
        f.writelines(['minimize_ligand                                              yes\n'])
        f.writelines(['simplex_max_iterations                                       1000\n'])
        f.writelines(['simplex_tors_premin_iterations                               0\n'])
        f.writelines(['simplex_max_cycles                                           1\n'])
        f.writelines(['simplex_score_converge                                       0.1\n'])
        f.writelines(['simplex_cycle_converge                                       1.0\n'])
        f.writelines(['simplex_trans_step                                           1.0\n'])
        f.writelines(['simplex_rot_step                                             0.1\n'])
        f.writelines(['simplex_tors_step                                            10.0\n'])
        f.writelines(['simplex_random_seed                                          0\n'])
        f.writelines(['simplex_restraint_min                                        no\n'])
        f.writelines(['atom_model                                                   all\n'])
        f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(config['dock6_path'])])
        f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(config['dock6_path'])])
        f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(config['dock6_path'])])
        f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
        f.writelines(['write_orientations                                           no\n'])
        f.writelines(['num_scored_conformers                                        1\n'])
        f.writelines(['rank_ligands                                                 no\n'])
    
    os.system('chmod 0700 {}'.format(run_sh_script))

    cmd = ['./{}'.format(run_sh_script)]
            
    return cmd

def scoring_finish_dock6_continuous_score(task, ret):
    try:
        with open('{}/ligand_out_scored.mol2'.format(task['tmp_run_dir']), 'r') as f:
            lines = f.readlines()
        
        task['score'] = float([x for x in lines[2].split(' ') if x][-1])
        task['status'] = "success"
    except:
        logger.error("failed parsing")    

    return

def scoring_start_mm_gbsa(task):
    config = load_config(config_path=task['config_path'])
    task['output_path'] = format_ligand(ligand_path=task['output_path'], file_format='mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    with open(run_sh_script, 'w') as f:
        # Getting Ligand Parameters: 
        f.writelines('export Chimera={}\n'.format(config['chimera_path']))
        f.writelines('charge=`$Chimera/bin/chimera --nogui --silent {} {}/charges.py`\n'.format(task['output_path'], task['tools_path']))
        f.writelines('antechamber -i {} -fi mol2 -o ligand_bcc.mol2 -fo mol2 -at gaff2 -c gas -rn LIG -nc $charge -pf y\n'.format(task['output_path']))
        f.writelines('parmchk2 -i ligand_bcc.mol2 -f mol2 -o ligand.frcmod\n')

        # Building Topology Files:
        f.writelines('tleap -f {}/tleap_r.in\n'.format(task['tools_path']))
        f.writelines('tleap -f {}/tleap_c.in\n'.format(task['tools_path']))
        
        # Run MD: 
        f.writelines('sander -O -i {}/min.in -p complex.prmtop -c complex.inpcrd -r min.rst -ref complex.inpcrd -o minim.out\n'.format(task['tools_path']))
        
        # Running MMPBSA.py
        f.writelines('MMPBSA.py -O -i {}/gbsa.in -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y  min.rst\n'.format(task['tools_path']))

    os.system(f'chmod 0700 {run_sh_script}')

    cmd = [f'./{run_sh_script}']

    return cmd

def scoring_finish_mm_gbsa(task, ret):
    try: 
        with open('{}/FINAL_RESULTS_MMPBSA.dat'.format(task['tmp_run_dir']), 'r') as f: 
            lines = f.readlines()

        lines = [x for x in lines if 'DELTA TOTAL' in x][0]
        
        task['score'] = float([x for x in lines.split(' ') if x][2])
        task['status'] = "success"
    except: 
        logger.error("failed parsing")
    
    return

def scoring_start_hawkins_gbsa(task):
    config = load_config(config_path=task['config_path'])
    task['output_path'] = format_ligand(ligand_path=task['output_path'], file_format='mol2')

    run_sh_script = os.path.join(task['tmp_run_dir'], "run.sh")
    with open(run_sh_script, 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(config['chimera_path'])])
        f.writelines(['export DOCK6={}\n'.format(config['dock6_path'])])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(config['receptor'])])
        f.writelines(['$DOCK6/bin/sphgen INSPH\n']) 
        f.writelines(['$DOCK6/bin/sphere_selector rec.sph {} 12.0 \n'.format(task['output_path'])])
        f.writelines(['$DOCK6/bin/showbox < box.in\n'])
        f.writelines(['$DOCK6/bin/grid -i grid.in\n'])
        f.writelines(['cd nchemgrid_GB\n'])
        f.writelines(['$DOCK6/bin/nchemgrid_GB\n'])
        f.writelines(['cd ../nchemgrid_SA\n'])
        f.writelines(['$DOCK6/bin/nchemgrid_SA\n'])
        f.writelines(['$DOCK6/bin/dock6 -i Hawkins_GBSA_Score.in\n'])

    os.system(f'chmod 0700 {run_sh_script}')  

    # Create INSPH File: 
    insph_path = os.path.join(task['tmp_run_dir'], "INSPH")
    with open(insph_path, 'w') as f: 
        f.writelines('rec.ms\n')
        f.writelines('R\n')
        f.writelines('X\n')
        f.writelines('0.0\n')
        f.writelines('4.0\n')
        f.writelines('1.4\n')
        f.writelines('rec.sph\n')

    # Create box.in File: 
    box_path = os.path.join(task['tmp_run_dir'], "box.in")
    with open(box_path, 'w') as f: 
        f.writelines('N\n')
        f.writelines('U\n')
        f.writelines('{}   {}    {}\n'.format(config['center_x'], config['center_y'], config['center_z']))
        f.writelines('{} {} {}\n'.format(config['size_x'], config['size_y'], config['size_z']))
        f.writelines('rec_box.pdb\n')

    grid_path = os.path.join(task['tmp_run_dir'], "grid.in")
    with open(grid_path, 'w') as f: 
        f.writelines('compute_grids                  yes\n')
        f.writelines('grid_spacing                   0.3\n')
        f.writelines('output_molecule                no\n')
        f.writelines('contact_score                  no\n')
        f.writelines('energy_score                   yes\n')
        f.writelines('energy_cutoff_distance         9999\n')
        f.writelines('atom_model                     a\n')
        f.writelines('attractive_exponent            6\n')
        f.writelines('repulsive_exponent             12\n')
        f.writelines('distance_dielectric            no\n')
        f.writelines('dielectric_factor              1\n')
        f.writelines('bump_filter                    yes\n')
        f.writelines('bump_overlap                   0.75\n')
        f.writelines('receptor_file                  {}\n'.format(config['receptor']))
        f.writelines('box_file                       rec_box.pdb\n')
        f.writelines('vdw_definition_file            {}/parameters/vdw_AMBER_parm99.defn\n'.format(config['dock6_path']))
        f.writelines('score_grid_prefix              solvent_grid\n')

    hawkins_gbsa_score_path = os.path.join(task['tmp_run_dir'], "Hawkins_GBSA_Score.in")
    with open(hawkins_gbsa_score_path, 'w') as f: 
        f.writelines('conformer_search_type                                        rigid\n')
        f.writelines('use_internal_energy                                          no\n')
        f.writelines('ligand_atom_file                                             {}\n'.format(task['output_path']))
        f.writelines('limit_max_ligands                                            no\n')
        f.writelines('skip_molecule                                                no\n')
        f.writelines('read_mol_solvation                                           no\n')
        f.writelines('calculate_rmsd                                               no\n')
        f.writelines('use_database_filter                                          no\n')
        f.writelines('orient_ligand                                                no\n')
        f.writelines('bump_filter                                                  no\n')
        f.writelines('score_molecules                                              yes\n')
        f.writelines('contact_score_primary                                        no\n')
        f.writelines('contact_score_secondary                                      no\n')
        f.writelines('grid_score_primary                                           no\n')
        f.writelines('grid_score_secondary                                         no\n')
        f.writelines('multigrid_score_primary                                      no\n')
        f.writelines('multigrid_score_secondary                                    no\n')
        f.writelines('dock3.5_score_primary                                        no\n')
        f.writelines('dock3.5_score_secondary                                      no\n')
        f.writelines('continuous_score_primary                                     no\n')
        f.writelines('continuous_score_secondary                                   no\n')
        f.writelines('footprint_similarity_score_primary                           no\n')
        f.writelines('footprint_similarity_score_secondary                         no\n')
        f.writelines('pharmacophore_score_primary                                  no\n')
        f.writelines('pharmacophore_score_secondary                                no\n')
        f.writelines('descriptor_score_primary                                     no\n')
        f.writelines('descriptor_score_secondary                                   no\n')
        f.writelines('gbsa_zou_score_primary                                       no\n')
        f.writelines('gbsa_zou_score_secondary                                     no\n')
        f.writelines('gbsa_hawkins_score_primary                                   yes\n')
        f.writelines('gbsa_hawkins_score_secondary                                 no\n')
        f.writelines('gbsa_hawkins_score_rec_filename                              rec_charged.mol2')
        f.writelines('gbsa_hawkins_score_solvent_dielectric                        78.5\n')
        f.writelines('gbsa_hawkins_use_salt_screen                                 no\n')
        f.writelines('gbsa_hawkins_score_gb_offset                                 0.09\n')
        f.writelines('gbsa_hawkins_score_cont_vdw_and_es                           no\n')
        f.writelines('gbsa_hawkins_score_grid_prefix                               solvent_grid\n')
        f.writelines('SASA_score_secondary                                         no\n')
        f.writelines('amber_score_secondary                                        no\n')
        f.writelines('minimize_ligand                                              no\n')
        f.writelines('atom_model                                                   all\n')
        f.writelines('vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(config['dock6_path']))
        f.writelines('flex_defn_file                                               {}/parameters/flex.defn'.format(config['dock6_path']))
        f.writelines('flex_drive_file                                              {}/parameters/flex_drive.tbl'.format(config['dock6_path']))
        f.writelines('ligand_outfile_prefix                                        gbsa_hawkins\n')
        f.writelines('write_orientations                                           no\n')
        f.writelines('num_scored_conformers                                        1\n')
        f.writelines('rank_ligands                                                 no\n')

    cmd = [f'./{run_sh_script}']

    return cmd

def scoring_finish_hawkins_gbsa(task, ret):
    with open(f"{task['tmp_run_dir']}/gbsa_hawkins_scored.mol2", 'r') as f: 
        lines = f.readlines()

    task['score'] = float([x for x in lines[2].split(' ') if x][-1])
    task['status'] = "success"
    
    return

DOCKING_PROGRAMS = {
    'MpSDockZN': {
        'start': docking_start_MpSDockZN,
        'end': docking_finish_MpSDockZN,
        'ligands': 'single',
    },
    'HDock': {
        'start': docking_start_HDock,
        'end': docking_finish_HDock,
        'ligands': 'single',
    },
    'dock6': {
        'start': docking_start_dock6,
        'end': docking_finish_dock6,
        'ligands': 'single',
    },
    'Flexx': {
        'start': docking_start_flexx,
        'end': docking_finish_flexx,
        'ligands': 'single',
    },
    'CovDock': {
        'start': docking_start_covdock,
        'end': docking_finish_covdock,
        'ligands': 'single',
    },
    'Glide_SP': {
        'start': docking_start_glide_sp,
        'end': docking_finish_glide_sp,
        'ligands': 'single',
    },
    'Glide_XP': {
        'start': docking_start_glide_xp,
        'end': docking_finish_glide_xp,
        'ligands': 'single',
    },
    'Glide_HTVS': {
        'start': docking_start_glide_htvs,
        'end': docking_finish_glide_htvs,
        'ligands': 'single',
    },
    'qvina02': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'qvina_w': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'vina': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'vina_carb': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'vina_xb': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'gwovina': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'AutodockVina_1.2': {
        'start': docking_start_vina,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'AutodockZN': {
        'start': docking_start_autodockzn,
        'end': docking_finish_vina,
        'ligands': "single"
    },
    'smina': {
        'start': docking_start_smina,
        'end': docking_finish_smina,
        'ligands': "single"
    },
    'adfr': {
        'start': docking_start_adfr,
        'end': docking_finish_adfr,
        'ligands': "single"
    },
    'plants': {
        'start': docking_start_plants,
        'end': docking_finish_plants,
        'ligands': "single"
    },
    'gnina': {
        'start': docking_start_gnina,
        'end': docking_finish_smina,
        'ligands': "single"
    },
    'rDock': {
        'start': docking_start_rdock,
        'end': docking_finish_rdock,
        'ligands': "single"
    },
    'M-Dock': {
        'start': docking_start_mdock,
        'end': docking_finish_mdock,
        'ligands': "single"
    },
    'MCDock': {
        'start': docking_start_mcdock,
        'end': docking_finish_mcdock,
        'ligands': "single"
    },
    'LigandFit': {
        'start': docking_start_ligandfit,
        'end': docking_finish_ligandfit,
        'ligands': "single"
    },
    'ledock': {
        'start': docking_start_ledock,
        'end': docking_finish_ledock,
        'ligands': "single"
    },
    'gold': {
        'start': docking_start_gold,
        'end': docking_finish_gold,
        'ligands': "single"
    },
    'iGemDock': {
        'start': docking_start_igemdock,
        'end': docking_finish_igemdock,
        'ligands': "single"
    },
    'idock': {
        'start': docking_start_idock,
        'end': docking_finish_idock,
        'ligands': "single"
    },
    'GalaxyDock3': {
        'start': docking_start_galaxydock3,
        'end': docking_finish_galaxydock3,
        'ligands': "single"
    },
    'autodock_cpu': {
        'start': docking_start_autodock_cpu,
        'end': docking_finish_autodock,
        'ligands': "single"
    },
    'autodock_gpu': {
        'start': docking_start_autodock_gpu,
        'end': docking_finish_autodock,
        'ligands': "single"
    },
    'autodock_koto': {
        'start': docking_start_autodock_koto,
        'end': docking_finish_autodock_koto,
        'ligands': "single"
    },
    'RLDock': {
        'start': docking_start_rldock,
        'end': docking_finish_rldock,
        'ligands': "single"
    },
    'PSOVina': {
        'start': docking_start_PSOVina,
        'end': docking_finish_PSOVina,
        'ligands': "single"
    },
    'LightDock': {
        'start': docking_start_LightDock,
        'end': docking_finish_LightDock,
        'ligands': "single"
    },
    'FitDock': {
        'start': docking_start_FitDock,
        'end': docking_finish_FitDock,
        'ligands': "single"
    },
    'Molegro': {
        'start': docking_start_Molegro,
        'end': docking_finish_Molegro,
        'ligands': "single"
    },
    'rosetta-ligand': {
        'start': docking_start_rosetta_ligand,
        'end': docking_finish_rosetta_ligand,
        'ligands': "single"
    },
    'SEED': {
        'start': docking_start_SEED,
        'end': docking_finish_SEED,
        'ligands': "single"
    },
    'FRED': {
        'start': docking_start_fred,
        'end': docking_finish_fred,
        'ligands': "single",
    },
    'scoring_nnscore2.0': {
        'start': scoring_start_nnscore2,
        'end': scoring_finish_nnscore2,
        'ligands': 'single',
    },
    'scoring_asp': {
        'start': scoring_start_gold_asp,
        'end': scoring_finish_gold_asp,
        'ligands': 'single',
    },
    'scoring_chemscore': {
        'start': scoring_start_gold_chemscore,
        'end': scoring_finish_gold_chemscore,
        'ligands': 'single',
    },
    'scoring_goldscore': {
        'start': scoring_start_gold_goldscore,
        'end': scoring_finish_gold_goldscore,
        'ligands': 'single',
    },
    'scoring_plp': {
        'start': scoring_start_gold_plp,
        'end': scoring_finish_gold_plp,
        'ligands': 'single',
    },
    'scoring_rf-score-vs': {
        'start': scoring_start_rf,
        'end': scoring_finish_rf,
        'ligands': 'single',
    },
    'scoring_smina': {
        'start': scoring_start_smina,
        'end': scoring_finish_smina,
        'ligands': 'single',
    },
    'scoring_start_gnina': {
        'start': scoring_start_gnina,
        'end': scoring_finish_gnina,
        'ligands': 'single',
    },
    'scoring_ad4': {
        'start': scoring_start_ad4,
        'end': scoring_finish_ad4,
        'ligands': 'single',
    },
    'scoring_vinardo': {
        'start': scoring_start_vinardo,
        'end': scoring_finish_vinardo,
        'ligands': 'single',
    },
    'scoring_vina': {
        'start': scoring_start_vina,
        'end': scoring_finish_vina,
        'ligands': 'single',
    },
    'scoring_PLANTS_chemplp': {
        'start': scoring_start_PLANTS_chemplp,
        'end': scoring_finish_PLANTS_chemplp,
        'ligands': 'single',
    },
    'scoring_PLANTS_plp': {
        'start': scoring_start_PLANTS_plp,
        'end': scoring_finish_PLANTS_plp,
        'ligands': 'single',
    },
    'scoring_PLANTS_plp95': {
        'start': scoring_start_PLANTS_plp95,
        'end': scoring_finish_PLANTS_plp95,
        'ligands': 'single',
    },
    'scoring_dock6': {
        'start': scoring_start_dock6_contact_score,
        'end': scoring_finish_dock6_contact_score,
        'ligands': 'single',
    },
    'scoring_dock6_continuous': {
        'start': scoring_start_dock6_continuous_score,
        'end': scoring_finish_dock6_continuous_score,
        'ligands': 'single',
    },
    'scoring_MM_GBSA': {
        'start': scoring_start_mm_gbsa,
        'end': scoring_finish_mm_gbsa,
        'ligands': 'single',
    },
    'scoring_Hawkins_GBSA': {
        'start': scoring_start_hawkins_gbsa,
        'end': scoring_finish_hawkins_gbsa,
        'ligands': 'single',
    },
}





def get_workunit_information():

    workunit_id = os.getenv('VFVS_WORKUNIT','') 
    subjob_id = os.getenv('VFVS_WORKUNIT_SUBJOB','')

    if(workunit_id == "" or subjob_id == ""):
        raise RuntimeError(f"Invalid VFVS_WORKUNIT and/or VFVS_WORKUNIT_SUBJOB")

    return workunit_id, subjob_id


def setup_job_storage_mode(ctx):

    ctx['job_storage_mode'] = os.getenv('VFVS_JOB_STORAGE_MODE', 'INVALID')

    if(ctx['job_storage_mode'] == "s3"):

        botoconfig = Config(
           region_name = os.getenv('VFVS_AWS_REGION'),
           retries = {
              'max_attempts': 50,
              'mode': 'standard'
           }
        )

        ctx['job_object'] = os.getenv('VFVS_CONFIG_JOB_OBJECT')
        ctx['job_bucket'] = os.getenv('VFVS_CONFIG_JOB_BUCKET')

        # Get the config information
        ctx['s3'] = boto3.client('s3', config=botoconfig)
    
    elif(ctx['job_storage_mode'] == "sharedfs"):
        ctx['job_tar'] = os.getenv('VFVS_CONFIG_JOB_TGZ')
    else:
        raise RuntimeError(f"Invalid jobstoragemode of {ctx['job_storage_mode']}. VFVS_JOB_STORAGE_MODE must be 's3' or 'sharedfs' ")


def get_subjob_config(ctx, workunit_id, subjob_id):

    if(ctx['job_storage_mode'] == "s3"):
        get_workunit_from_s3(ctx, workunit_id, subjob_id, 
            ctx['job_bucket'], ctx['job_object'], ctx['temp_dir'])
    elif(ctx['job_storage_mode'] == "sharedfs"):
        get_workunit_from_sharedfs(ctx, workunit_id, subjob_id,
            ctx['job_tar'], ctx['temp_dir'])
    else:
        raise RuntimeError(f"Invalid jobstoragemode of {ctx['job_storage_mode']}. VFVS_JOB_STORAGE_MODE must be 's3' or 'sharedfs' ")




def process(ctx):


    ctx['vcpus_to_use'] = int(os.getenv('VFVS_VCPUS', 1))
    ctx['run_sequential'] = int(os.getenv('VFVS_RUN_SEQUENTIAL', 0))

    # What job are we running?

    workunit_id, subjob_id =  get_workunit_information()

    # Setup paths appropriately depending on if we are using S3
    # or a shared FS
    
    setup_job_storage_mode(ctx)

    # This includes all of the configuration information we need
    # After this point ctx['main_config'] has the configuration options
    # and we have specific subjob information in ctx['subjob_config']

    get_subjob_config(ctx, workunit_id, subjob_id)

    # Update some of the path information

    process_config(ctx)

    ctx['workunit_id'] = workunit_id
    ctx['subjob_id'] = subjob_id

    ctx.pop('s3', None)

    print(ctx['temp_dir'])

    ligand_format = ctx['main_config']['ligand_library_format']


    # Need to expand out all of the collections in this subjob
    subjob = ctx['subjob_config']

    metadata = {
        'workunit_id': workunit_id,
        'subunit_id': subjob_id,
        'ligand_library_format': ligand_format,
        'vcpus_to_use': ctx['vcpus_to_use']
    }


    download_queue = Queue()
    unpack_queue = Queue()
    collection_queue = Queue()
    docking_queue = Queue()
    summary_queue = Queue()
    upload_queue = Queue()



    try:
        downloader_processes = []
        for i in range(0, math.ceil(ctx['vcpus_to_use'] / 8.0)):
            downloader_processes.append(Process(target=downloader, args=(download_queue, unpack_queue, summary_queue, ctx['temp_dir'])))
            downloader_processes[i].start()

        unpacker_processes = []
        for i in range(0, math.ceil(ctx['vcpus_to_use'] / 8.0)):
            unpacker_processes.append(Process(target=untar, args=(unpack_queue, collection_queue)))
            unpacker_processes[i].start()

        collection_processes = []
        for i in range(0, math.ceil(ctx['vcpus_to_use'] / 8.0)):
            # collection_process(ctx, collection_queue, docking_queue, summary_queue)
            collection_processes.append(Process(target=collection_process, args=(ctx, collection_queue, docking_queue, summary_queue)))
            collection_processes[i].start()

        docking_processes = []
        for i in range(0, ctx['vcpus_to_use']):
            # docking_process(docking_queue, summary_queue)
            docking_processes.append(Process(target=docking_process, args=(ctx, docking_queue, summary_queue)))
            docking_processes[i].start()

        # There should never be more than one summary process
        summary_processes = []
        summary_processes.append(Process(target=summary_process, args=(ctx, summary_queue, upload_queue, metadata)))
        summary_processes[0].start()

        uploader_processes = []
        for i in range(0, 2):
            # docking_process(docking_queue, summary_queue)
            uploader_processes.append(Process(target=upload_process, args=(ctx, upload_queue)))
            uploader_processes[i].start()


        for collection_key in subjob['collections']:
            collection = subjob['collections'][collection_key]
            collection_name, collection_number = collection_key.split("_", maxsplit=1)
            collection['collection_number'] = collection_number
            collection['collection_name'] = collection_name

            # Checking if dynamic tranche filtering is enabled and prescreen_mode is disabled (since the collections names in the prescreen mode do not have tranche names with properties encoded)
            if (int(ctx['main_config']['prescreen_mode']) == 0) and (int(ctx['main_config']['dynamic_tranche_filtering']) == 1):

                # Compiling the user's regex
                try:
                    compiled_regex = re.compile(ctx['main_config']['dynamic_tranche_filtering_regex'])
                except re.error as e:
                    print(f"Invalid regular expression: {e}")
                    return None

                match = re.search(compiled_regex, collection_name)

                if not match:
                    print(f" Skipping collection {collection_name} due to dynamic tranche filtering. Tranche of collection: {collection_name}, regex:{ctx['main_config']['dynamic_tranche_filtering_regex']}")
                    #logging.error(f"Tranche {collection_name} not contained in part of the library specified by regular expression. Skipping.")
                    break

            download_item = {
                'collection_key': collection_key,
                'collection': collection,
                'ext': "tar.gz",
            }

            download_queue.put(download_item)

            # Don't overflow the queues
            while download_queue.qsize() > 25:
                time.sleep(0.2)

        flush_queue(download_queue, downloader_processes, "download")
        flush_queue(unpack_queue, unpacker_processes, "unpack")
        flush_queue(collection_queue, collection_processes, "collection")
        flush_queue(docking_queue, docking_processes, "docking")
        flush_queue(summary_queue, summary_processes, "summary")
        flush_queue(upload_queue, uploader_processes, "upload")
    except Exception as e:
        logging.error(f"Received exception {e}, terminating")

        for process in [*downloader_processes, *unpacker_processes, *collection_processes, *docking_processes, *summary_processes, *uploader_processes]:
            print(process.name)
            process.kill()
        sys.exit(1)



def flush_queue(queue, processes, description):
    logging.error(f"Sending {description} flush")
    for process in processes:
        queue.put(None)
    logging.error(f"Join {description}")
    for process in processes:
        process.join()
        if(process.exitcode != 0):
            raise RuntimeError(f'Process from {description} exited with {process.exitcode}')
    logging.error(f"Finished Join of {description}")



def main():

    ctx = {}

    log_level = os.environ.get('VFVS_LOGLEVEL', 'INFO').upper()
    logging.basicConfig(level=log_level)

    ctx['tools_path'] = os.getenv('VFVS_TOOLS_PATH', "/opt/vf/tools/bin")

    # Temp directory information
    temp_path = os.getenv('VFVS_TMP_PATH', None)
    if(temp_path):
        temp_path = os.path.join(temp_path, '')

    with tempfile.TemporaryDirectory(prefix=temp_path) as temp_dir:
        ctx['temp_dir'] = temp_dir

        # stat = shutil.disk_usage(path)
        stat = shutil.disk_usage(ctx['temp_dir'])
        if(stat.free < (1024 * 1024 * 1024 * 1)):
            raise RuntimeError(f"VFVS needs at least 1GB of space free in tmp dir ({ctx['temp_dir']}) free: {stat.free} bytes")


        print(ctx['temp_dir'])
        process(ctx)
        print(ctx['temp_dir'])


if __name__ == '__main__':
    main()
