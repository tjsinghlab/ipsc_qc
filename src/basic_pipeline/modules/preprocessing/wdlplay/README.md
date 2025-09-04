# wdlplay

## Overview
This is a Python package that can submit, monitor, and manage WDL workflows on a Caper server or Cromwell server on-prem (`cromwell-singh01`). 
## Main Features

## Installation
Clone this repository. In the directory of the repo locally, run `make create-env` to get a clean environment. Then run `make install` to install the wdlplay package and requirements. 

## Usage
We have three options to run Caper/Cromwell server:
1) Run `srun --cpus-per-task=8 --mem=32G --pty bash` to get an interactive node to start the Caper server on. 
2) You can use the `login-singh` node alternatively, but this is not recommended since it is a shared node. 

If you use these two methods, in the node, run `screen` to start a screen session (alternatively, use `tmux` if you are more familiar with that). 

3) You can use the already set up `cromwell-singh01` on cluster. This does not require any set up, just use `cromwell-singh01` as the server hostname later on. 

In the screen session, do `module load caper singularity`.

### Start mysql server (optional, when you need to run hundreds on thousands of jobs)

Pick a port to run your mysql server on, and add the following to your `~/.caper/default.conf`:
```
db=mysql
mysql-db-port=$MYSQL_PORT
```
Then, run
```
bash src/sandbox/run_mysql_server.sh /gpfs/commons/groups/singh_lab/users/<user>/mysql-server $MYSQL_PORT
```
The argument `/gpfs/commons/groups/singh_lab/users/<user>/mysql-server` needs to be an empty folder that mysql server is going to write to.

### Start Caper server
Choose another port for the Caper server. It's recommended to choose unique port numbers (e.g. 8345) to avoid port clashes. Then run:
```
caper server \
--port $CAPER_PORT \
--mysql-db-port $MYSQL_PORT \
--backend-file /gpfs/commons/groups/singh_lab/resources/caper/gcp-backend.conf
```
After seeing `Ready to take submissions`, the Caper server is up and running.

### Submit Jobs Using wdlplay
Make a subfolder in `src` for your project. Use `src/plink-template` as a reference.

Put both your WDL script and a `run.py` (or called `sandbox.py`) in `src/PROJECT`. Run `run.py` to submit jobs. In this script, initialize wdlplayer first, passing in the hostname and port arguments. `hostname` should be the node you started Caper server on (or `cromwell-singh01` if you used that), and `port` should be the port of the Caper server (8000 for `cromwell-singh01`). 
```
runner = wdlplayer(
    outdir=outdir,
    filedir=filedir,
    localdata=localdir,
    tmpdir=tmpdir
    caper_backend_tag='slurm',
    hostname='pe2-cc2-066',
    port=8345)
```
Specify the name of the job, and path to your WDL script. Prepare inputs to the WDL workflow in a dictionary. 
```
name = 'plink-template'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/repos/wdlplay/src/plink-template/plink-template.wdl'

input_wdl = {
    "plink_workflow.input_dir": '/gpfs/commons/groups/singh_lab/users/dongwang/plink_test/',
    "plink_workflow.input_prefix": 'gpc_AFR_bip_i_plink',
    "plink_workflow.docker": '/gpfs/commons/groups/singh_lab/resources/images/genetics_latest.sif'
}
```
Then run the following to submit the job. 
```
outfile = runner.submit(name=name, wdl=wdl, inputs=input_wdl, server=True, no_deepcopy=True)
```

### Debugging Workflows
If a workflow fails, run the following to get the error of the workflow:
```
runner.debug('eff9ced8-3b9e-4d04-95d0-6230150fca57') # pass in a workflow id
```

## WDL Script Design
Refer to `src/plink-template/plink-template.wdl` for an example WDL script.
### Runtime Attributes
```
cpus: 8
cpu: 8
mem: "16G"
memory: "16G"
singularity: docker
```
Specify both `cpus` and `cpu`, both `mem` and `memory` to accommodate for different backends. Having all options there will not lead to errors.

For the Docker/Singularity image, however, you can only specify one type of image. 
* If your backend in NYGC cluster, Singularity is recommended. Use `singularity: path_to_singularity_image`. The path should be a local path to a `.sif` file to optimize environment initialization speed.
* If your backend is GCP VM instance, only Docker is allowed. Use `docker: path_to_docker_image`. The Docker path should be a `gcr.io` path pointing to a Docker image in Google Cloud Registry. 

## Cite

## Maintainer

singhlab@nygenome.org

## Acknowledgements

## Release Notes