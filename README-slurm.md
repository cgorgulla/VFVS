

## Getting Started with VirtualFlow Virtual Screening (VFVS) with Slurm


### Prerequisites

VFVS requires Python 3.x (it has been tested with Python 3.9.4). Additionally, it requires packages to be installed that are not included at many sites. As a result, we recommend that you create a virtualenv for VFVS that can be used to run jobs.

Make sure that virtualenv is installed:
```bash
python3 -m pip install --user --upgrade virtualenv
```

Create a virtualenv (this example creates it under $HOME/vfvs_env):
```bash
python3 -m virtualenv $HOME/vfvs_env
```

Enter the virtualenv and install needed packages
```bash
source $HOME/vfvs_env/bin/activate
python3 -m pip install boto3 pandas pyarrow jinja2
```

If the `vfvs_get_top_results.py` script shall be used to extract the top results after screening, additional python packages are needed. To install these:
```bash
python3 -m pip install tqdm
```

To exit a virtual environment:
```bash
deactivate
```

When running VFVS commands, you should enter the virtualenv that you have setup using:
```bash
source $HOME/vfvs_env/bin/activate
```

As noted later, you will want to include this virtualenv as part of the Slurm script that VFVS runs for a job.


### Install VFVS

#### Download the VFVS Code

Login to a node that is part of your cluster (generally the login node) and execute the following to obtain the latest version of the code. This should be placed in a location that is on a shared filesystem that is available on all of the compute nodes.

```bash
git clone https://github.com/VirtualFlow/VFVS.git
cd VFVS
```

#### Update the configuration file

The file is in `tools/templates/all.ctrl` and the options are documented in the file itself. A few key Slurm-based parameters:

Job Configuration:

- `batchsystem`: Set this to `slurm` if you are running with the Slurm Workload Manager scheduler
- `threads_per_docking`: This is how many threads should be run per docking execution. This is almost always '1' since VFVS will run multiple docking executions in parallel for higher efficiency vs more threads per single docking.
- `threads_to_use`: Set this to the number of threads cores that a single job should use. Many Slurm clusters are configured so a single job consumes an entire node; if that is the case, set this value to the number of cores on the compute node that VFVS will be run on. (This will generally be the same as `slurm_cpus` below.)
- `program_timeout`: Seconds to wait until deciding that a single docking execution has timed out.
- `job_storage_mode`: When using Slurm, this must be set to `sharedfs`. This setting means that the data (input collections and output) is on a shared filesystem that can be seen across all nodes. This may be a scratch directory.

Slurm-specific Configuration:

- `slurm_template`: This is the path to the template file that will be used to submit the invidual jobs that VFVS submits to slurm. Use this to set specific options that may be required by your site. This could include setting an account number, specific partition, etc. Note: Please update this to include the virtualenv that you setup earlier (if needed)
- `slurm_array_job_size`: How many jobs can be run in a single array job. (The default for a Slurm scheduler is 1000, so unless your site has changed it then the limit will likely be this)
- `slurm_array_job_throttle`: VFVS uses Slurm array jobs to more efficiently submit the jobs.  This setting limits how many jobs from a single job array run could be run at the same time. Setting `slurm_array_job_throttle` to `slurm_array_job_size` means there will be no throttling.
- `slurm_partition`: Name of slurm partition to submit to
- `slurm_cpus`: Number of compute cores that should be requested per job.

Job-sizing:

- `dockings_todo_per_queue`: This determines how many ligands should be processed at a minimum per job. A value of '10000' would mean that each subjob with `slurm_cpus` number of CPUs should dock this number of ligands prior to completing. In general jobs should run for approximately 30 minutes or more. How long each docking takes depends on the receptor, ligand being docked, and docking program-specific settings (such as `exhaustiveness`). Submitting a small job to determine how long a docking will take is often a good idea to size these before large runs.


The ligands to be processed should be included in the file within `tools/templates/todo.all`. This file can be automatically generated from the VirtualFlow website.


### Data for Virtual Screening

The location of the collection files (containing the list of ready-to-dock ligands) to be used in the screening should be located in `collection_folder` (defined in `all.ctrl`). 
Apart from the input ligands, all other inputs used for docking (such as the receptor file and docking config file) are located as subdirectories within the `docking_scenario_basefolder`.

VFVS expects that collection data will be stored in one of two different directory structures, defined as ``hash`` or ``metatranche``. Typically this will be the `metatranche` setting.


#### `metatranche` addressing

In this setting the collections are stored in a format where the collection string “`AECCAACEACDB`” would be stored under the prefix: `AECC/AACE/ACDB-0000`. This data format is the default for the datasets currently available on the VirtualFlow website. Note, if you have downloaded data from the VirtualFlow website, the format it will download in will require extraction of data to expand the first level of tar archives.


#### `hash` addressing

This evenly distributes files across different prefixes, which can be beneficial for various filesystems and metadata. This is most appropriate for situations where there may be 1B+ ligands to process. With this addressing mode, the collection string “`AECCAACEACDB`” would be stored under `65/0b/<datasetname>/<datatype>/AECCAACEACDB/0000000.tar.gz`. The hash address (`650b`) is generated from the first 4 characters of a stringified SHA256 hexdigest of `AECCAACEACDB/0000000`.

```bash
[ec2-user ~]$ echo -n "AECCAACEACDB/0000000" | sha256sum | cut -c1-4
650b
```

This output format can be automatically generated from newer versions of VFLP. (If you already have a dataset in the `metatranche` mode it is not recommended to transform it into the `hash` addressing mode.)

### Run a Job

Job commands should be run directly from within the `tools` directory of the VFVS package.

#### Prepare Workflow

```bash
cd tools
./vfvs_prepare_folders.py
```

If you have previously setup a job in this directory the command will let you know that it already exists. If you are sure you want to delete the existing data, then run with `--overwrite`.

Once you run this command the workflow is defined using the state of `all.ctrl` and `todo.all` at that time. Changes to those files at this point will not be used unless `vfvs_prepare_folders.py` is run again.

#### Generate Workunits

VFVS can process billions of ligands, and in order to process these efficiently it is helpful to segment this work into smaller chunks. A workunit is a segment of work that contains many 'subjobs' that are the actual execution elements. Often a workunit will have many subjobs and each subjob will contain about 60 minutes worth of computation.

```bash
./vfvs_prepare_workunits.py
```

Pay attention to how many workunits are generated. The final line of output will provide the number of workunits generated.

#### Submit the job to run on Slurm

The following command will submit workunits 1 and 2. The default configuration with
Slurm will use 200 subjobs per workunit, so this will submit 2x200 (400) subjobs
(assuming that each workunit was full). Each subjob takes `slurm_cpus` cores when running.

How long each job takes will be dependent on the parameters that were set as part of the `all.ctrl` and the docking scenarios themselves.

```bash
./vfvs_submit_jobs.py 1 2
```

Once submitted, the jobs will be visible in the Slurm queue (`squeue`)

#### Monitor Progress

The progress of the jobs can be monitored using the Slurm scheduler commands. 
For debugging, the slurm errors are contained within `../workflow/workunits`.


## Viewing Results

The results of the job will be placed in the location specified within the `all.ctrl`.
By default, this is saved in: `/output-files/{workunit-number}`. 

#### `metatranche` addressing

Assuming the `metatranche` addressing is being used, the data can be found in the bucket specified under the `workload` directory as part of the VFVS installation.


````
workload/output/<scenario_name>/results/<metatranche>/<tranche>/<collection_num>.tar.gz
workload/output/<scenario_name>/logfiles/<metatranche>/<tranche>/<collection_num>.tar.gz
workload/output/<scenario_name>/summaries/<metatranche>/<tranche>/<collection_num>.tar.gz
workload/output/ligand-lists/<metatranche>/<tranche>/<collection_num>.status.gz
workload/output/ligand-lists/<metatranche>/<tranche>/<collection_num>.json.gz
````

The top results can also be extracted using the `tools/vfvs_get_top_results.py` script (execute `vfvs_get_top_results.py --help` for more info).

