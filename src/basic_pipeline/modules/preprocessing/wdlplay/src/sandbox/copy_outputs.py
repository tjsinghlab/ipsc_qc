import os
import json
from logzero import logger

import wdlplay
from wdlplay.wdlplayer import wdlplayer

import argparse

def main():
    # parser = argparse.ArgumentParser()
    
    # parser.add_argument('--start', type=int, help='Start parameter')
    # parser.add_argument('--end', type=int, help='End parameter')
    
    # args = parser.parse_args()
    
    # start_param = args.start
    # end_param = args.end

    # print(start_param, end_param)

    runner = wdlplayer(
        outdir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay-logs',
        filedir='/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay-logs',
        caper_backend_tag='slurm',
        current_run=1)

    metadata_dir = "/gpfs/commons/groups/singh_lab/users/dongwang/wdlplay/src/sandbox/log/0.1.0"
    metadata_list = sorted(os.listdir(metadata_dir))
    metadata_list = [os.path.join(metadata_dir,m) for m in metadata_list]

    runner.copy_outputs("/gpfs/commons/groups/singh_lab/projects/hbccseq/data/processed-WGS", metadata_list=metadata_list)

if __name__ == "__main__":
    main()