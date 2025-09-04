# %%
#region
import os
import json
from logzero import logger

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath

from datatracker import *
from wdlplay.wdlplayer import wdlplayer

#endregion

os.environ['VERSION'] = '0.1.0'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

#%%
tr = Tracker()
entry = Entry(tag='warp-parallel',
              description='testing submit parallel warp pipeline',
              category='Processing',
              module='Processing')

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
    current_run=1)
    # cromwell_backend_file='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/backend.conf')
#endregion


# %%
#---------------------------------------------------
name = 'warp-non-cmc'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq-project/warp-working/WholeGenomeReprocessing_Pipeline.wdl'
# inputs = runner.get_inputs(wdl)

input_files_path = "/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/WGSprocessing_input_files_non_CMC.json"
with open(input_files_path, 'r') as json_file:
    input_files = json.load(json_file)

input_settings_path = "/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/WGSprocessing_input_settings.json"
with open(input_settings_path, 'r') as json_file:
    input_settings = json.load(json_file)

inputs = [{**i, **input_settings} for i in input_files]
print(inputs)
print(len(inputs))
names = [name + "-" + str(i+1) for i in range(len(inputs))]

# %%
outfiles = [runner.submit(name=n, wdl=wdl, inputs=i, server=True, port=8200) for n,i in zip(names, inputs)]

# %%
runner.update_metadata()

# %%
runner.resubmit(wdl=wdl, server=True, port=8200, keyword="warp", keyword_exclude="re")
