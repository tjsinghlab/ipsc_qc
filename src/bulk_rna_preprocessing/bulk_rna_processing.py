# %%

#wdl script sourced from warp-pipelines/rnaseq-processing/wdl (sdwang008)

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
entry = Entry(tag='eqtl',
              description='prep eQTL',
              category='Analysis',
              module='eQTL')

# %%
is_cloud=False
set_cloud(cloud=is_cloud)

#region
runner = wdlplayer(
    outdir='/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025_2',
    filedir=filedir,
    localdata=localdir,
    tmpdir='/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025_2/.caper_tmp/',
    caper_backend_tag='slurm',
    hostname='login-singh',
    port=8200,
    gcp_project='singh-comp-d-271c')
#endregion

#%%
name = 'process_bulk_rna'
wdl = '/gpfs/commons/groups/singh_lab/users/dongwang/hbccseq/src/RNA-seq/1.0_process-raw/RNAseq_pipeline_fastq.wdl'

#%%
directory = "/gpfs/commons/groups/singh_lab/projects/bd2village/data/hart_data_may2025"
onlyfiles = os.listdir(directory)
onlyfiles = [f for f in onlyfiles if ".fastq" in f or ".fq" in f]
onlyfiles = sorted(onlyfiles)

inputs_list = []
for i in range(0, len(onlyfiles), 2):
    inputs_list.append((
        [os.path.join(directory, onlyfiles[i]), os.path.join(directory, onlyfiles[i+1])], # first element of the tuple, which is a list of the R1 and R2 fastq files
        onlyfiles[i].split("_")[0] # second element of the tuple, which is the prefix of the fastq files
        # CHANGE HOW THE PREFIX IS EXTRACTED BASED ON YOUR FILE NAMING SCHEME
        ))
# for loop through the list of fastq files at an increment of 2. For each pair of fastq files, make a tuple where the first element is a list of the R1 and R2 fastq files and the second element is the prefix of the fastq files

inputs = [{
      "rnaseq_pipeline_fastq_workflow.fastqs": i[0], # list of length 2, the first element is the R1 fastq and the second element is the R2 fastq
      "rnaseq_pipeline_fastq_workflow.prefix": i[1], # sample name of the fastq files
      "rnaseq_pipeline_fastq_workflow.star_index_oh75": "/gpfs/commons/groups/singh_lab/resources/RNAseq/star_index_oh75", # reference files for STAR
      "rnaseq_pipeline_fastq_workflow.rsem_reference": "/gpfs/commons/groups/singh_lab/resources/RNAseq/rsem_reference", # reference files for RSEM
      "rnaseq_pipeline_fastq_workflow.genes_gtf": "/gpfs/commons/groups/singh_lab/resources/RNAseq/gencode.v39.GRCh38.genes.collapsed_only.gtf", # reference GTF
      "rnaseq_pipeline_fastq_workflow.outdir": "/gpfs/commons/groups/singh_lab/users/kjakubiak/ipsc/ron_data_may2025_2/" # where to save outputs to
    } for i in inputs_list]
names = [name + "-" + str(i+1) for i in range(len(inputs))]

#%%
outfile = runner.submit(name=names[0], wdl=wdl, inputs=inputs[0], server=True, no_deepcopy=False)
# outfile = [runner.submit(name=n, wdl=wdl, inputs=i, server=True, no_deepcopy=False) for n,i in zip(names, inputs)]
# %%

outfiles = []
failed_samples = []

# Define the function to check if the sample has already been processed
def is_sample_processed(output_dir, sample_name):
    # Example condition: check for an expected output file
    output_file = os.path.join(output_dir, f"{sample_name}.processed")  # Modify as needed
    return os.path.exists(output_file)

# Initial loop to submit workflows
for name, input_set in zip(names, inputs):
    output_dir = "/gpfs/commons/groups/singh_lab/users/kjakubiak/hbcc_snRNAseq/pool_analysis_2"
    sample_name = input_set["cellbender_remove_background.run_cellbender_remove_background_gpu.sample_name"]
    
    if is_sample_processed(output_dir, sample_name):
        logger.info(f"Skipping {name}: output already exists.")
        continue  # Skip this sample
    
    try:
        logger.info(f"Submitting workflow for: {name}")
        outfile = runner.submit(
            name=name,
            wdl=wdl,
            inputs=input_set,
            server=True,
            no_deepcopy=False
        )
        outfiles.append(outfile)
    except Exception as e:
        logger.error(f"Error submitting workflow for {name}: {e}")
        failed_samples.append((name, input_set))

# Retry loop for failed submissions
while failed_samples:
    logger.info(f"Retrying {len(failed_samples)} failed samples...")
    remaining_failed = []

    for name, input_set in failed_samples:
        try:
            logger.info(f"Retrying workflow for: {name}")
            outfile = runner.submit(
                name=name,
                wdl=wdl,
                inputs=input_set,
                server=True,
                no_deepcopy=False
            )
            outfiles.append(outfile)
        except Exception as e:
            logger.error(f"Error resubmitting workflow for {name}: {e}")
            remaining_failed.append((name, input_set))
    
    # Update the list of failed samples
    failed_samples = remaining_failed

logger.info(f"All workflows processed. Total successful: {len(outfiles)}")

# %%


