# %%
import os
import json
from logzero import logger
import argparse

import wdlplay
from findhere import here, init_directories, set_cloud, reldir, relpath
from wdlplay.wdlplayer import wdlplayer

# ------------------------------------------------
# CLI arguments
# ------------------------------------------------
parser = argparse.ArgumentParser(description="Run WDL stage 2: Germline variant calling")
parser.add_argument("--output_dir", required=True, help="Output directory")
parser.add_argument("--sample", required=True, help="Sample name")
args = parser.parse_args()

# ------------------------------------------------
# Versioning & logging
# ------------------------------------------------
os.environ['VERSION'] = '0.1.0'
log_dir = os.path.join(args.output_dir, "logs", os.environ['VERSION'], args.sample)
os.makedirs(log_dir, exist_ok=True)
os.environ['WDLPLAY_LOGDIR'] = log_dir

# ------------------------------------------------
# Metadata & directories
# ------------------------------------------------
cloudir, localdir, filedir = init_directories(__file__)
is_cloud = False
set_cloud(cloud=is_cloud)

# ------------------------------------------------
# Runner setup
# ------------------------------------------------
runner = wdlplayer(
    outdir=args.output_dir,
    filedir=filedir,
    localdata=localdir,
    tmpdir=os.path.join(args.output_dir, "wdlplay_logs", "WholeGenomeProcessing", ".caper_tmp"),
    caper_backend_tag='local'
)

# ------------------------------------------------
# Inputs
# ------------------------------------------------
name = f"rna-variant-calling-gvcf-{args.sample}"
#wdl = "gatk4-rna-best-practices.wdl"
wdl = os.path.join(filedir, 'gatk4-rna-best-practices.wdl')


inputs_template = "/pipeline/modules/preprocessing/wdlplay/warp-pipelines/GATK_variant_calling/gatk4-rna-germline-variant-calling.inputs.json"
inputs = json.load(open(inputs_template))

# BAM + BAI from MarkDuplicates
bam_dir = os.path.join(args.output_dir, "Mark_duplicates_outputs")
bam_file = os.path.join(bam_dir, f"{args.sample}.Aligned.sortedByCoord.out.md.bam")
bai_file = bam_file + ".bai"

if not os.path.exists(bam_file) or not os.path.exists(bai_file):
    raise FileNotFoundError(f"Missing BAM or BAI for sample {args.sample} in {bam_dir}")

# Variant calling output dir
vcf_out = os.path.join(args.output_dir, "variant_calling", args.sample)
os.makedirs(vcf_out, exist_ok=True)

inputs["RNAseq.inputBam"] = bam_file
inputs["RNAseq.inputBamIndex"] = bai_file
inputs["RNAseq.outdir"] = vcf_out
inputs["RNAseq.sampleNameCustom"] = args.sample

# ------------------------------------------------
# Submit workflow
# ------------------------------------------------
logger.info(f"Submitting WDL workflow for sample {args.sample}")
outfile = runner.submit(
    name=name,
    wdl=wdl,
    inputs=inputs,
    server=False,
    no_deepcopy=False,
    check_metadata=False
)

logger.info(f"Workflow submitted. Output: {outfile}")
# %%