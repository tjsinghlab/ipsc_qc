# %%
#region
import os
import json
from logzero import logger

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath

from datatracker import *
from wdlplay.wdlplayer import wdlplayer

import pandas as pd
#endregion

os.environ['VERSION'] = '0.1.0'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

#%%
tr = Tracker()
entry = Entry(tag='template',
              description='plink template',
              category='Processing',
              module='Processing')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir=os.getcwd(),
    filedir=filedir,
    localdata=localdir,
    tmpdir=os.getcwd(),
    caper_backend_tag='slurm',
    gcp_project='singh-comp-d-271c',
    hostname='login-singh',
    port=8000)
#endregion


# %%
#---------------------------------------------------
name = 'plink-template'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/repos/wdlplay/src/plink-template/plink-template.wdl'

input_wdl = {
    "plink_workflow.input_dir": '/gpfs/commons/groups/singh_lab/users/dongwang/plink_test/',
    "plink_workflow.input_prefix": 'gpc_AFR_bip_i_plink',
    "plink_workflow.docker": '/gpfs/commons/groups/singh_lab/resources/images/genetics_latest.sif'
}

# %%

outfiles = runner.submit(name=name, wdl=wdl, inputs=input_wdl, server=True)

# %%
runner.debug('8f4217ad-794f-4ada-a1b1-717f947159f7')