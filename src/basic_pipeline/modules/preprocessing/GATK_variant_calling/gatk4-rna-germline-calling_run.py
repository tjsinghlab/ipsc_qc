# %%

#code sourced from warp-pipelines/rna-germline-variant-calling (sdwang008)
#region
import os
import json
from logzero import logger

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath

from datatracker import *
from wdlplay.wdlplayer import wdlplayer

from google.cloud import storage

#endregion

os.environ['VERSION'] = '0.1.0'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

#%%
tr = Tracker()
entry = Entry(tag='rna-germline-variant-calling',
              description='Call germline variants from RNA-seq data using GATK4 best practices.',
              category='processing',
              module='processing')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir='output_dir/wdlplay_outputs/WholeGenomeProcessing/',
    filedir=filedir,
    localdata=localdir,
    tmpdir='output_dir/wdlplay_logs/WholeGenomeProcessing/.caper_tmp/',
    caper_backend_tag='slurm',
    hostname='',
    port=8200,
    gcp_project='')
#endregion

#%%
name = 'rna-variant-calling-gvcf'
wdl = 'gatk4-rna-best-practices.wdl'

#%%
# construct inputs list
inputs = json.load(open('gatk4-rna-germline-variant-calling.inputs.json'))

bam_dir = "/path/to/star/bams/Mark_duplicates_outputs"
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
outfile = [runner.submit(name=n, wdl=wdl, inputs=i, server=True, no_deepcopy=False, check_metadata=False) for n,i in zip(names, inputs_list)]
# %%