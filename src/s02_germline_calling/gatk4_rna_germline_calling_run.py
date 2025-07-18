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
entry = Entry(tag='sortBAM',
              description='sort BAM',
              category='processing',
              module='processing')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay_logs/WholeGenomeProcessing',
    filedir=filedir,
    localdata=localdir,
    tmpdir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay_logs/WholeGenomeProcessing/.caper_tmp/',
    caper_backend_tag='slurm',
    hostname='login-singh',
    port=8200,
    gcp_project='singh-comp-d-271c')
#endregion

#%%
name = 'rna-variant-calling'
wdl = 'gatk4-rna-best-practices.wdl'

#%%
# map sample names
fastqs = os.listdir("/gpfs/commons/groups/singh_lab/projects/bd2village/data/BD2_Dec_2024")
fastqs = [f for f in fastqs if "R1_001.fastq.gz" in f]
mapping = {f.split("_")[1]: f.split("_")[0] for f in fastqs}

#%%
# construct inputs list
inputs = json.load(open('gatk4-rna-germline-variant-calling.inputs.json'))

bam_dir = "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/bd2_village/wdl_pipeline/outs/Mark_duplicates_outputs/"
bams = sorted(os.listdir(bam_dir))
bams = [b for b in bams if b.endswith(".bam")]

inputs_list = []
for i in range(0, len(bams)):
    current_inputs = inputs.copy()
    current_inputs["RNAseq.inputBam"] = os.path.join(bam_dir, bams[i])
    current_inputs["RNAseq.inputBamIndex"] = os.path.join(bam_dir, bams[i]) + ".bai"
    current_inputs["RNAseq.outdir"] = "/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/RNA-variant-calling/"
    current_inputs["RNAseq.sampleNameCustom"] = mapping[bams[i].split(".")[0]]
    inputs_list.append(current_inputs)

names = [name + "-" + str(i+1) for i in range(len(inputs_list))]

#%%
outfile = [runner.submit(name=n, wdl=wdl, inputs=i, server=True, no_deepcopy=False) for n,i in zip(names, inputs_list)]
# %%
