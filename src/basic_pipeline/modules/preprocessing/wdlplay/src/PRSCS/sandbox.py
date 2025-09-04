# %%
#region
import os
import json
from logzero import logger

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath

from datatracker import *
from wdlplay.wdlplayer import wdlplayer
import hail as hl
from joblib import Parallel, delayed

import pandas as pd
#endregion

os.environ['VERSION'] = '0.1.0'
cloudir, localdir, filedir = init_directories(__file__)#, relbase='src')

#%%
tr = Tracker()
entry = Entry(tag='prscs',
              description='prscs calculation',
              category='Processing',
              module='Processing')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay_logs/PRSCS',
    filedir=filedir,
    localdata=localdir,
    tmpdir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay_logs/PRSCS/.caper_tmp/',
    caper_backend_tag='slurm',
    hostname='login-singh',
    port=8200,
    gcp_project='singh-comp-d-271c')
    # cromwell_backend_file='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/config/backend.conf')
#endregion


# %%
#---------------------------------------------------
name = 'PRSCS'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/repos/wdlplay/src/PRSCS/prscs.wdl'

# %%
input_prscs = {
    "PRSCS.raw_sumstats": '',
    "PRSCS.tag": '',
    "PRSCS.bim": '/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/data/WGS/5.0_prs/1.0_hbcc_prs_NO_HWE.bim',
    "PRSCS.sample_size": '',
    "PRSCS.outdir": '/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/data/WGS/5.0_prs/2.0_hbcc_prscs_grch37/'
}

metadata_df = pd.read_table('gs://singh-comp-d-271c-scoreregistry/data/processing/3.0_generate-munge/0.1.58/registry-munged.tsv', sep='\t')

tags = ['2018_adhd', '2019_mdd', '2020_scz', '2020_bip', '2019_an', '2019_asd', '2018_ocd', '2018_neuroticism', '2020_cud', '2018_scz_bd', '2022_suicide_attempt', '2021_bip_all', '2023_mdd', '2023_adhd', '2022_scz', '2023_suicide_eur', '2021_bip_noUKBB', '2021_bip_bdi', '2021_bip_bdii', '2023_suicide_eur_no_ukb']
# tags = ['2020_bip','2019_mdd']
metadata_df_subset = metadata_df[metadata_df['tag'].isin(tags)]


# %%

sumstat_num = 1
for index, row in metadata_df_subset.iterrows():
    if row['exclude_from_analysis'] == "TRUE":
        continue
    tag = row['tag']
    N = row['N']
    raw_sumstats=row['processed_tsv_path']
    if pd.isna(raw_sumstats):
        continue

    input_prscs['PRSCS.tag'] = tag
    input_prscs['PRSCS.sample_size'] = N     
    input_prscs['PRSCS.raw_sumstats'] = raw_sumstats

    outfile = runner.submit(name=name+str(sumstat_num), wdl=wdl, inputs=input_prscs, server=True, no_deepcopy=False)

    sumstat_num += 1
            
    # concat_job = concat_and_append_header_job(b, infiles=chrom_jobs, header_list=['CHR', 'SNP', 'BP', 'A1', 'A2', 'BETA'])
    # concat_job = gzip_chain_job(concat_job)
    # b.write_output(concat_job.ofile, outfile)

# %%
# runner.update_metadata()

# # %%
# runner.resubmit(wdl=wdl, server=True, port=8200)

# # %%
# runner.resume()

# # %%
# runner.get_error()
