#find wdl script at https://github.com/tjsinghlab/hbccseq/blob/main/src/analysis/2.0_eQTL/2.0_peer_factors.wdl
# %%
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
entry = Entry(tag='peer',
              description='peer factors',
              category='Analysis',
              module='eQTL')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir='/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/peer',
    filedir=filedir,
    localdata=localdir,
    tmpdir='/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/peer/.caper_tmp/',
    caper_backend_tag='slurm',
    hostname='ne1dc4-006',
    port=8200,
    gcp_project='singh-comp-d-271c')
#endregion

#%%
name = 'peer_factors'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/src/analysis/2.0_eQTL/2.0_peer_factors.wdl'
region = 'iPSC'
inputs = {"peer_factors_workflow.phenotype_file": "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/lognorm_for_peer.tsv",
          "peer_factors_workflow.prefix": region,
          "peer_factors_workflow.num_peer": 30,
          "peer_factors_workflow.outdir": "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/peer/"
          }

#%%
outfile = runner.submit(name=name+"_"+region, wdl=wdl, inputs=inputs, server=True, no_deepcopy=False)
# %%
