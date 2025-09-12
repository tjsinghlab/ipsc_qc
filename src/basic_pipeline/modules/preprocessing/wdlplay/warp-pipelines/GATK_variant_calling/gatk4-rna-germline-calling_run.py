# %%

#code sourced from warp-pipelines/rna-germline-variant-calling (sdwang008)
#region
import os
import json
from logzero import logger

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath
import argparse
from datatracker import *
from wdlplay.wdlplayer import wdlplayer

from google.cloud import storage

#endregion

# ------------------------------------------------
# CLI arguments
# ------------------------------------------------
parser = argparse.ArgumentParser(description="Run WDL stage 2: Germline variant calling")
parser.add_argument("--output_dir", required=True, help="Top-level output directory (writable, host-mounted)")
parser.add_argument("--sample", required=True, help="Sample name")
args = parser.parse_args()

# ------------------------------------------------
# Versioning & logging
# ------------------------------------------------
os.environ['VERSION'] = '0.1.0'

# Logs live inside the user-provided output_dir
log_dir = os.path.join(args.output_dir, "logs", os.environ['VERSION'], args.sample)
os.makedirs(log_dir, exist_ok=True)
os.environ['WDLPLAY_LOGDIR'] = log_dir

# ------------------------------------------------
# Metadata & directories
# ------------------------------------------------
cloudir, localdir, filedir = init_directories(__file__)

# Optional tracking (disabled for now)
# tr = Tracker()
# entry = Entry(
#     tag='rna-germline-variant-calling',
#     description='Call germline variants from RNA-seq data using GATK4 best practices.',
#     category='processing',
#     module='processing'
# )

is_cloud = False
set_cloud(cloud=is_cloud)

# ------------------------------------------------
# Runner
# ------------------------------------------------
runner = wdlplayer(
    outdir=os.path.join(args.output_dir, "wdlplay_outputs", "WholeGenomeProcessing"),
    filedir=filedir,
    localdata=localdir,
    tmpdir=os.path.join(args.output_dir, "wdlplay_logs", "WholeGenomeProcessing", ".caper_tmp"),
    caper_backend_tag='slurm'
    # hostname='', port=8200, gcp_project=''  # add if needed
)

#%%
name = 'rna-variant-calling-gvcf'
wdl = 'gatk4-rna-best-practices.wdl'

#%%
# construct inputs list
inputs = json.load(open('gatk4-rna-germline-variant-calling.inputs.json'))

bam_dir = "${outdir}Mark_duplicates_outputs"
bams = sorted(os.listdir(bam_dir))
bams = [b for b in bams if b.endswith(".bam")]

inputs_list = []
for i in range(0, len(bams)):
    current_inputs = inputs.copy()
    current_inputs["RNAseq.inputBam"] = os.path.join(bam_dir, bams[i])
    current_inputs["RNAseq.inputBamIndex"] = os.path.join(bam_dir, bams[i]) + ".bai"
    current_inputs["RNAseq.outdir"] = "/output_dir/variant_calling/" # make this empty directory first
    current_inputs["RNAseq.sampleNameCustom"] = bams[i].split(".")[0]
    inputs_list.append(current_inputs)
inputs_list = [inputs_list[i] for i in [1]]

names = [name + "-" + str(i+1) for i in range(len(inputs_list))]

#%%
outfile = [runner.submit(name=n, wdl=wdl, inputs=i, server=False, no_deepcopy=False, check_metadata=False) for n,i in zip(names, inputs_list)]
# %%