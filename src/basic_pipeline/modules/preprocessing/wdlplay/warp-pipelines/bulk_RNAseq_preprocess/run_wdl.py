#!/usr/bin/env python3
import os
import json
import argparse
from logzero import logger
from wdlplay.wdlplayer import wdlplayer
from findhere import init_directories, set_cloud
from datatracker import Tracker, Entry

# ------------------------------------------------
# CLI arguments
# ------------------------------------------------
parser = argparse.ArgumentParser(description="Run WDL stage 1: RNA-seq preprocessing")
parser.add_argument("--fastq1", required=True, help="Path to FASTQ R1 file")
parser.add_argument("--fastq2", required=True, help="Path to FASTQ R2 file")
parser.add_argument("--output_dir", required=True, help="Output directory")
parser.add_argument("--sample", required=True, help="Sample name (prefix)")
args = parser.parse_args()

os.makedirs(args.output_dir, exist_ok=True)

# ------------------------------------------------
# Metadata & tracking (kept from original)
# ------------------------------------------------
os.environ['VERSION'] = '0.1.0'
cloudir, localdir, filedir = init_directories(__file__)

tr = Tracker()
entry = Entry(
    tag='eqtl',
    description='prep eQTL',
    category='Analysis',
    module='eQTL'
)

is_cloud = False
set_cloud(cloud=is_cloud)

runner = wdlplayer(
    outdir=args.output_dir,
    filedir=filedir,
    localdata=localdir,
    tmpdir=os.path.join(args.output_dir, ".caper_tmp"),
    caper_backend_tag='slurm',
    hostname='',
    port=8200,
    gcp_project=''
)

# ------------------------------------------------
# Define workflow + inputs
# ------------------------------------------------
name = f"process_bulk_rna-{args.sample}"
wdl = "/pipeline/modules/preprocessing/wdlplay/warp-pipelines/bulk_RNAseq_preprocess/RNAseq_pipeline_fastq.wdl"  # path inside container

inputs = {
    "rnaseq_pipeline_fastq_workflow.fastqs": [args.fastq1, args.fastq2],
    "rnaseq_pipeline_fastq_workflow.prefix": args.sample,
    "rnaseq_pipeline_fastq_workflow.star_index_oh75": "/ref/star_index_oh75",
    "rnaseq_pipeline_fastq_workflow.rsem_reference": "/ref/rsem_reference",
    "rnaseq_pipeline_fastq_workflow.genes_gtf": "/ref/gencode.v39.GRCh38.genes.collapsed_only.gtf",
    "rnaseq_pipeline_fastq_workflow.outdir": args.output_dir
}

# ------------------------------------------------
# Submit workflow
# ------------------------------------------------
logger.info(f"Submitting WDL workflow for sample {args.sample}")
outfile = runner.submit(
    name=name,
    wdl=wdl,
    inputs=inputs,
    server=False,
    no_deepcopy=False
)

logger.info(f"Workflow submitted. Output: {outfile}")
