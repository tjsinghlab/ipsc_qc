"""
This script submits WDL jobs to cromwell server on SLURM to run LDSC on processed
sumstats from scoreregistry. The results are then parsed and saved to a tsv file.
"""
# %%
#region
import os
import json
from logzero import logger

from findhere import *

from datatracker import *
from wdlplay.wdlplayer import wdlplayer
import hail as hl
from joblib import Parallel, delayed

import pandas as pd
import numpy as np
#endregion

os.environ['VERSION'] = '0.1.3'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

#%%
# tr = Tracker()

# entry = Entry(tag='mtag-prscs',
#         description='Run prscs calculation on mtag outage sumstats',
#         category='Analysis',
#         module='MTAG',
#         version=os.environ['VERSION']) 
# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
outdir=f"/gpfs/commons/groups/singh_lab/users/fshiau/project/h2enrich/data/ldscs_liability/{os.environ['VERSION']}/ldsc_h2_outdir"
tmpdir = "/gpfs/commons/groups/singh_lab/users/fshiau/repos/wdlplay/.caper_tmp/.caper_server/"
runner = wdlplayer(
    outdir=os.path.join(outdir, "wdlplay_logs"),
    filedir=filedir,
    localdata=outdir,
    tmpdir=tmpdir,
    caper_backend_tag='slurm',
    current_run=1,
    gcp_project='singh-comp-d-271c')#,
#    config='/gpfs/commons/groups/singh_lab/resources/caper/gcp-backend.conf')
    # cromwell_backend_file='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/backend.conf')
#endregion
#%%
os.makedirs(runner.outdir, exist_ok=True)
# %%
#---------------------------------------------------
name = 'LDSC'
wdl = '/gpfs/commons/groups/singh_lab/users/fshiau/repos/wdlplay/src/ldsc/ldsc.wdl'
#%%

metadata = pd.read_csv(
    "gs://singh-comp-d-271c-scoreregistry/data/processing/2.0_standardize-sumstats/0.1.54/registry-processed.tsv",
    sep="\t")
#%%
atlas = pd.read_csv("gs://singh-comp-d-271c-scoreregistry/data/metadata/metadata.csv",
                    sep = ",")
#%%
metadata.loc[:,"popprev"]=metadata.apply(
    lambda x: atlas.loc[atlas.Trait.str.lower().str.contains(
        x.Trait.lower()),"popprev"].to_list() if x.Trait==x.Trait else None,
        axis=1)
#%%
for i, row in enumerate(metadata.popprev):
    prev = [p for p in row if p==p]
    prev = np.average(prev) if len(prev) > 0 else None
    metadata.loc[i, "popprev"] = prev
# %%
# runner.project = 'singh-comp-d-271c'

#%%
for i, row in metadata.iterrows():
    if row["exclude_from_analysis"] == True:
        continue
    elif os.path.exists(os.path.join(outdir, f"{row['tag']}_h2.log")):
        metadata.loc[i, "h2_log_liability"] = os.path.join(outdir, f"{row['tag']}_h2.log")
        continue
    elif pd.isnull(row.popprev):
        continue
    elif np.isnan(row.popprev):
        continue
    else:
        input_ldsc = {
            "LDSC.raw_sumstats_gz": row["processed_tsv_path"],
            "LDSC.merge_alleles": "/gpfs/commons/groups/singh_lab/projects/h2enrich/data/resources/w_hm3.snplist",
            "LDSC.sample_size": row["N"],
            "LDSC.ldsc_h2_outdir": outdir,
            "LDSC.phenotype": row["tag"],
            "LDSC.ancestry" : "EUR",
            "LDSC.panel": "1KG",
            "LDSC.liability": "liability",
            "LDSC.samprev": f"{round(row['Ncase']/row['N'], 3):.3f}",
            "LDSC.poprev": str(row["popprev"])}
        outfile = runner.submit(name=f"ldsc_{i}", wdl=wdl,
                                inputs=input_ldsc, server=True, port=8205)
        metadata.loc[i, "h2_log_liability"] = os.path.join(outdir, f"{row['tag']}_h2.log")
    

#%%
metadata.to_csv(os.path.join(outdir, "registry-processed.tsv"), sep="\t", index=False)
# %%
from h2enrich.logparser import logparser_h2
# %%
h2_parse = []
for i, row in metadata.iterrows():
    if row["exclude_from_analysis"] == True:
        continue
    elif pd.isnull(row["popprev"]):
        continue
    else:
        h2_parse += [logparser_h2(row["tag"],  row["h2_log_liability"])]
# %%
h2 = pd.concat(h2_parse).reset_index(drop = "index")
#%%
h2.loc[:,"popprev"] = metadata.set_index("tag").loc[h2.tag,"popprev"].values
# %%
h2.to_csv(os.path.join(outdir, "h2_liability_results.tsv"), sep="\t", index=False)
# %%
h2.to_csv(f"gs://singh-h2enrich/data/ldscs/{os.environ['VERSION']}/h2_liability_results.tsv",
          sep="\t", index=False)
# %%
