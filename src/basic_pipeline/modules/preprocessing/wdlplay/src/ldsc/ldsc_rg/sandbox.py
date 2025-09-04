"""
This script submits WDL jobs to cromwell server on SLURM to run LDSC genetic 
correlation between processed sumstats from scoreregistry. The results are 
then parsed and saved to a tsv file.

**NOTE**
When starting caper server, set `--max-concurrent-tasks` to int greater than
n_phenotypes^2-1 to avoid overloading the server.
"""
# %%
# region
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
# endregion

os.environ['VERSION'] = '0.1.5'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

# %%
# tr = Tracker()

# entry = Entry(tag='mtag-prscs',
#         description='Run prscs calculation on mtag outage sumstats',
#         category='Analysis',
#         module='MTAG',
#         version=os.environ['VERSION'])
# %%
is_cloud=False
set_cloud(cloud=is_cloud)

# region
outdir=f"/gpfs/commons/groups/singh_lab/users/fshiau/project/h2enrich/data/"
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
# endregion
# %%
os.makedirs(runner.outdir, exist_ok=True)
# %%
# ---------------------------------------------------
name = 'LDSC'
wdl = '/gpfs/commons/groups/singh_lab/users/fshiau/repos/wdlplay/src/ldsc/ldsc_rg/ldsc_rg.wdl'
# %%

metadata = pd.read_csv(
    "gs://singh-comp-d-271c-scoreregistry/data/processing/2.0_standardize-sumstats/0.1.57/registry-processed.tsv",
    sep="\t")
# %%
metadata = pd.concat(
    [
        metadata,
        pd.DataFrame(
            {
                "gwastype": "gpc",
                "processed_tsv_path": [
                    "gs://gpc_array/data/analysis/10.0_hl_gwas/10.5.0_hl_process-bip-sumstats/0.2.43/gpc-AFR-bipi-processed_sumstat.tsv.gz",
                    "gs://gpc_array/data/analysis/10.0_hl_gwas/10.5.0_hl_process-bip-sumstats/0.2.43/gpc-AMR-bipi-processed_sumstat.tsv.gz",
                    "gs://gpc_array/data/analysis/10.0_hl_gwas/10.5.0_hl_process-bip-sumstats/0.2.43/gpc-EUR-bipi-processed_sumstat.tsv.gz",
                ],
                "N": ["4666", "3644", "5769"],
                "tag": ["gpc-AFR-bipi", "gpc-AMR-bipi", "gpc-EUR-bipi"],
            }
        ),
    ],
    ignore_index=True,
)
# %%
metadata = metadata.iloc[:5,:]
# %%
phenotypes = metadata.loc[metadata.exclude_from_analysis != True, "tag"].to_list()
# %%
input_ldsc = {
    "LDSC.raw_sumstats_gz": metadata.loc[
        metadata.exclude_from_analysis != True, "processed_tsv_path"
    ].to_list(),
    "LDSC.merge_alleles": "/gpfs/commons/groups/singh_lab/projects/h2enrich/data/resources/w_hm3.snplist",
    "LDSC.sample_size": metadata.loc[
        metadata.exclude_from_analysis != True, "N"
    ].to_list(),
    "LDSC.munge_outdir": os.path.join(outdir, "munged"),
    "LDSC.rg_outdir": os.path.join(outdir, "rg"),
    "LDSC.phenotypes": phenotypes,
    "LDSC.panel": "1KG",
}
outfile = runner.submit(
    name=f"test", wdl=wdl, inputs=input_ldsc, server=True, port=8205
)

# %%
metadata.to_csv(os.path.join(outdir, "registry-processed.tsv"), sep="\t", index=False)
# %%
from h2enrich.logparser import logparser_rg
# %%
h2_parse = []
for pheno1 in phenotypes:
    for pheno2 in phenotypes:
        try:
            df_tmp = logparser_rg(
                f"{pheno1}_{pheno2}",
                os.path.join(outdir, "rg",f"{pheno1}_{pheno2}.log")
                )
            df_tmp.loc[:,"tag1"] = pheno1
            df_tmp.loc[:,"tag2"] = pheno2
            h2_parse += [df_tmp]
        except:
            print(f"{pheno1}_{pheno2} failed")
# %%
h2 = pd.concat(h2_parse).reset_index(drop = "index")
# %%
h2.to_csv(os.path.join(outdir, "rg_results.tsv"), sep="\t", index=False)
# %%
h2.to_csv(f"gs://singh-h2enrich/data/ldscs/{os.environ['VERSION']}/rg_results.tsv",
          sep="\t", index=False)
# %%
