import json
import os
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor
import sys
import time
from collections import Counter
import pandas as pd

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class wdlplayer(object):
    def __init__(self, outdir, filedir, localdata=None, tmpdir=None, config=None, caper_backend_tag='local', cromwell_backend_file=None, inputs=None, gcp_project=None, gcp_zone=None, current_run=1, hostname=None, port=None):
        """
        Initializes a WDLPlayer object to run and submit wdl jobs.

        Parameters
        ----------
        outdir : str
            The path to the output directory. If used with the`findhere` package, 
            this should be set to `cloudir` for GCP runs or `localdir` for local runs.
        filedir : str
            The path to the file directory.
        localdata : str, optional
            The path to the local data directory. If used with the`findhere` package, 
            this should be set to `localdir`. Default is None.
        tmpdir : str, optional
            The path to the temporary directory. 
            Default is None, which sets tmpdir to `$outdir/.caper_tmp`
        config : str, optional
            The path to the Caper's backend config JSON file. Default is None, which uses caper defaults
        inputs : str, optional
            The path to the inputs file. Default is None.
        project : str, optional
            The name of the Google Cloud project. Default is None.
        zone : str, optional
            The name of the Google Cloud zone. Default is None.

        Examples
        --------
        runner = wdlplayer(
            outdir=local,
            filedir='/path/to/files',
            localdata='/path/to/local/data',
            tmpdir='/path/to/tmp',
            config='/path/to/config.json',            inputs='/path/to/inputs.txt',
            project='my_project',
            zone='my_zone'
        )
        """
        # Assigns outdir and tmpdir
        self.outdir = outdir
        self.filedir = filedir
        self.localdata = localdata if localdata is not None else self.outdir
        localdata = self.localdata
        self.tmpdir = tmpdir if tmpdir is not None else os.path.join(
            self.outdir, '.caper_tmp')
        
        # Define version
        self.version = 'debug' if 'VERSION' not in os.environ else os.environ['VERSION']
        
        logdir_env = os.environ.get("WDLPLAY_LOGDIR", None)
        if logdir_env:
            self.logdir = logdir_env
        else:
            self.logdir = os.path.join(self.outdir, "log", self.version)

        os.makedirs(self.logdir, exist_ok=True)
            
        # Assigns config json path
        self.caper_backend_tag = caper_backend_tag
        self.cromwell_backend_file = cromwell_backend_file
        self.config = config
        self.metadata = None
        self.uid = None
        
        # Assigns cloud parameters
        self.project = gcp_project
        self.zone = gcp_zone
        
        #---------------------------------------------------
        # Same as hailrunner
        # to keep syntax, change outdir
        #self.outdir = filedir

        

        # Define hostname and port
        self.hostname = hostname

        self.port = port

        # Create outdir and basic .jobs file structure
        if not os.path.isdir(self.outdir) or not os.path.isdir(os.path.join(self.outdir, 'log', self.version)):
            os.makedirs(os.path.join(self.outdir, "log", self.version) + '/')
        self.logdir = os.path.join(self.outdir, 'log', self.version)

        if not os.path.isdir(os.path.join(self.outdir, 'script', self.version)):
            os.makedirs(os.path.join(
                self.outdir, "script", self.version) + '/')
        self.scriptdir = os.path.join(self.outdir, 'script', self.version)
        
        # add current links
        if os.path.isdir(os.path.join(self.outdir, "current")):
            shutil.rmtree(os.path.join(self.outdir, "current"))
        os.makedirs(os.path.join(self.outdir, "current"))
        os.symlink(self.logdir, os.path.join(
            self.outdir, 'current', 'log'), target_is_directory=True)
        os.symlink(self.scriptdir, os.path.join(
            self.outdir, 'current', 'script'), target_is_directory=True)
        
        #add local /data folder links
        if localdata != '':
            if os.path.islink(os.path.join(self.outdir, 'data')):
                os.remove('data')
            os.symlink(localdata, os.path.join(
                self.outdir, 'data'), target_is_directory=True)

        # 0-based counter for submits made in script. Used for naming and tracking.
        self.current_run = current_run

        # 0-based counter that stores a particular run count for restoration later
        self._setr = []
        self._sig_popr = None

        # Initialize run-specific array attributes of PyRunner class
        self.current_job = []

        #---------------------------------------------------

        # Read in inputfofn can be a list of files or read from a file
        self.fofn = []
        if isinstance(inputs, str):
            self.fofn.append(inputs)
        else:
            self.fofn.append('')

        #---------------------------------------------------
        # fix outdir
        self.outdir = outdir

        self.get_jars()

    def get_jars(self):
        lines = [line.strip() for line in open(
            os.path.join(os.environ['HOME'], '.caper/default.conf'), 'r')]
        self.cromwell_jar = [line.split('=')[1]
                            for line in lines if re.search('^cromwell', line)][0]
        self.womtool_jar = [line.split('=')[1]
                            for line in lines if re.search('^womtool', line)][0]

    def submit(self, name, wdl, inputs, options=None, imports=None, server=False, cromwell_backend_file=None, caper_backend_tag=None, silent=False, dry_run=False, no_deepcopy=True, check_metadata=True):
    
        #---------------------------------------------------
        # hailrunner (modded from submit_job)
        if len(self.fofn) == self.current_run:
            self.fofn.append([])
        self.current_job = [wdl, inputs, options]
        self.fofn[self.current_run] = self.outdir
    
        #---------------------------------------------------
        # checkpointing
        log_path = os.path.join(
            self.logdir, "run{}.{}.metadata.json".format(self.current_run, name))
        
        if os.path.exists(log_path) and check_metadata:
            self.metadata = json.load(open(log_path, 'r'))
            self.uid = self.metadata['id']
            self.update_metadata(uid=self.uid, log_path=log_path)
            self.metadata = json.load(open(log_path, 'r'))

            if self.metadata['status'] == 'Succeeded':
                logger.info('Workflow {}.{} successfully completed in an earlier attempt.'.format(self.current_run, name))
            elif self.metadata['status'] == 'Running':
                logger.info('Workflow {}.{} is still running. Its cromwell id is {}.'.format(self.current_run, name, self.metadata['id']))
            elif self.metadata['status'] == 'Failed':
                if server:
                    logger.info('Workflow {}.{} failed. Its cromwell id is {}. Consider using the resubmit() function.'.format(self.current_run, name, self.metadata['id']))
                elif self.metadata['status'] != 'Succeeded':
                    raise Exception(f'Job {self.uid} has not succeeded. See log file at {log_path}')   
            return             
            
        else:            
            #---------------------------------------------------
            # Specify script and log path
            script_path = os.path.join(
                self.scriptdir,
                "run{}.{}".format(self.current_run, name)
            )
            
            if os.path.isfile(self.current_job[0]):
                shutil.copy2(self.current_job[0], f'{script_path}.wdl')
            if isinstance(inputs, str) and os.path.isfile(self.current_job[1]):
                shutil.copy2(self.current_job[1], f'{script_path}.inputs.json')
            if isinstance(options, str) and os.path.isfile(self.current_job[2]):
                shutil.copy2(self.current_job[2], f'{script_path}.options.json')
            
            if isinstance(inputs, dict):
                json.dump(inputs, open(
                    f'{script_path}.inputs.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=4)
                inputs = f'{script_path}.inputs.json'
                
            if isinstance(options, dict):
                json.dump(options, open(
                    f'{script_path}.options.json', 'w', encoding='utf-8'), ensure_ascii=False, indent=4)
                options = f'{script_path}.options.json'
                    
            #---------------------------------------------------
            # command
            cmd = ['caper']
            
            # SERVER or RUN mode
            if server:
                cmd.extend(['submit'])
            else:
                cmd.extend(['run'])

            # specify wdl, input JSON and label (name)
            cmd.extend([wdl,
                        '--inputs', inputs,  
                        '--str-label', name
                        ])

            # specify options JSON if not None
            if isinstance(options, str):
                cmd.extend(['--options', options])
                
            # specify imports.zip if not None
            if isinstance(imports, str):
                cmd.extend(['--imports', imports])

            # dry run command
            if dry_run:
                cmd.extend(['--dry-run'])
            
            # disable deepcopy
            if no_deepcopy:
                cmd.extend(['--no-deepcopy'])

            if isinstance(self.hostname, str):
                cmd.extend(['--hostname', self.hostname])

            # define caper backend tag
            caper_backend_tag = caper_backend_tag if caper_backend_tag is not None else self.caper_backend_tag

            # define Caper pipeline options using command or config file (default, not Cromwell)
            if self.config is not None:
                # use config to specify RUN or SERVER mode
                cmd.extend(['--conf', self.config])
            elif not server:
                # specify backend options in RUN mode
                # otherwise, need to be scpecified in PORT mode
                if caper_backend_tag == 'gcp':
                    cmd.extend([
                        '--backend', 'gcp',
                        '--gcp-prj', self.project,
                        '--gcp-zones', self.zone,
                        '--gcp-out-dir', self.outdir,
                        '--gcp-loc-dir', self.tmpdir,
                        '--use-google-cloud-life-sciences'
                    ])
                elif caper_backend_tag == 'slurm':
                    cmd.extend(['--backend', 'slurm'])
                    cmd.extend(['--local-out-dir', self.outdir])
                    cmd.extend(['--local-loc-dir', self.tmpdir])
                elif caper_backend_tag == 'local':
                    cmd.extend(['--backend', 'local'])
                    cmd.extend(['--local-out-dir', self.outdir])
                    cmd.extend(['--local-loc-dir', self.tmpdir])
                else:
                    raise ValueError("caper_backend_tag must be one of 'gcp', 'slurm', or 'local'.")
            elif server:
                # specify PORT to connect in SERVER mode
                # other backend options are defined in the server itself
                if isinstance(self.port, int):
                    cmd.extend(['--port', str(self.port)])
                if caper_backend_tag == 'gcp':
                    cmd.extend(['--gcp-loc-dir', self.tmpdir])
            
            # specify Cromwell config if not None
            if isinstance(cromwell_backend_file, str):
                cmd.extend(['--backend-file', cromwell_backend_file])
            elif self.cromwell_backend_file is not None:
                cmd.extend(['--backend-file', self.cromwell_backend_file])
            
            # cromwell log file
            if not server:
                cmd.extend(['--cromwell-stdout',
                            os.path.join(self.logdir, "run{}.{}.cromwell.out".format(self.current_run, name))])

            #---------------------------------------------------
            # execute
            output = self._execute(cmd, silent)
            
            if dry_run:
                return(output)
            
            if not server:
                outfile = self._parse_run_output(output)
                self._process_job_outputs(outfile, log_path)
            elif server:
                # return uid without "completing" job
                self.uid = self._parse_server_output(output)
                print('Pausing to let backend write metadata...')
                if check_metadata:
                    time.sleep(5)
                    self.metadata = self.get_metadata(self.uid)
                    self.write_json(log_path, self.metadata)
                return(self.metadata)
        
        #---------------------------------------------------
        # done (from hailrunner)
        self.done(outfile)
        
        return(outfile)

    def resume(self, keyword=None, keyword_exclude=None):
        '''resumes wdlplayer and runs update_metadata() and resubmit() in a while loop'''
        self.update_metadata(keyword=keyword,keyword_exclude=keyword_exclude)
        while self.failures(keyword=keyword,keyword_exclude=keyword_exclude):
            self.update_metadata()
            self.resubmit(keyword=keyword,keyword_exclude=keyword_exclude)
            time.sleep(600)

    def failures(self, keyword=None, keyword_exclude=None):
        job_ids, statuses, names = self.list(keyword=keyword,keyword_exclude=keyword_exclude)
        return "Failed" in statuses
    
    def update_metadata(self, uid=None, log_path=None, keyword=None, keyword_exclude=None):
        if uid is None and log_path is None:
            job_ids, statuses, names = self.list(keyword=keyword,keyword_exclude=keyword_exclude)
            for uid, name in zip(job_ids, names):
                log_path = os.path.join(self.logdir, "run{}.{}.metadata.json".format(self.current_run, name))
                metadata = self.get_metadata(uid)
                self.write_json(log_path, metadata)
        else:
            metadata = self.get_metadata(uid)
            name = metadata['labels']['caper-str-label']
            log_path = os.path.join(self.logdir, "run{}.{}.metadata.json".format(self.current_run, name))
            self.write_json(log_path, metadata)

    def resubmit(self, wdl, server=True, port=8200, keyword=None, keyword_exclude=None):
        '''check for statuses of current jobs on server and resubmit failed ones'''
        job_ids, statuses, names = self.list(keyword=keyword,keyword_exclude=keyword_exclude)

        failed = [i for i in range(len(statuses)) if statuses[i] == "Failed"]
        failed_jobs = []
        retry_names = []
        for f in failed:
            failed_jobs.append(job_ids[f])
            i = 1
            retry_name = names[f]+"-re"+str(i)
            log_path = os.path.join(self.logdir, "run{}.{}.metadata.json".format(self.current_run, retry_name))
            while os.path.exists(log_path):
                i += 1
                retry_name = names[f]+"-re"+str(i)
                log_path = os.path.join(self.logdir, "run{}.{}.metadata.json".format(self.current_run, retry_name))
            retry_names.append(retry_name) 

        for i in range(len(failed_jobs)):
            metadata = self.get_metadata(uid=failed_jobs[i])
            retry_inputs = json.loads(metadata["submittedFiles"]["inputs"])
            self.submit(name=retry_names[i], wdl=wdl, inputs=retry_inputs, server=server, port=port)

    def copy_outputs(self, copy_dst, fofn_directory, keyword=None, keyword_exclude=None, metadata_list=None, output_keys=None):
        if not os.path.exists(copy_dst):
            os.makedirs(copy_dst)

        # if list of metadata is not provided, get metadata from caper server
        if keyword is not None or keyword_exclude is not None:
            job_ids, statuses, names = self.list(keyword=keyword,keyword_exclude=keyword_exclude)
            metadata_list = [self.get_metadata(j) for j in job_ids]
        
        # with a list of metadata files, look through the metadata files and get output paths
        if metadata_list is not None:
            output_paths = {}
            fofn = []
            fofn_paths = []
            sample_names = []
            for m in metadata_list:
                metadata = json.load(open(m, "r"))
                sample_name = json.loads(metadata['submittedFiles']['inputs'])['WholeGenomeReprocessing.base_file_name']
                sample_names.append(sample_name)
                sample_dst = os.path.join(copy_dst,sample_name)
                if not os.path.exists(sample_dst):
                    os.makedirs(sample_dst)

                output_path = self.get_output_directory(metadata=metadata, sample_dst=sample_dst)

                sample_path = {}
                for k,v in output_path.items():
                    sample_path.update({k: list(v.values())[0]})
                    output_paths.update(v)
                fofn_paths.append(sample_path)
                print(sample_name)

            fofn_df = pd.DataFrame(fofn_paths)
            fofn_df.index = sample_names
            fofn_df.to_csv(fofn_directory)
            
            origins = list(output_paths.keys())
            dsts = list(output_paths.values())

            origin_paths_file = os.path.join(copy_dst,"origins.txt")
            dst_paths_file = os.path.join(copy_dst,"dsts.txt")
            with open(origin_paths_file, 'w') as origins_file:
                origins_file.write('\n'.join(origins))
            with open(dst_paths_file, 'w') as dsts_file:
                dsts_file.write('\n'.join(dsts))

            xargs_cmd = "paste -d \" \" {} {} | xargs -P 8 -n 2 cp".format(origin_paths_file, dst_paths_file)
            # print(xargs_cmd)
            subprocess.run(xargs_cmd, shell=True)

    def get_output_directory(self, metadata, sample_dst, keys=None):
        # job status check, then get output paths
        if metadata['status'] == 'Running':
            logger.info('This job is still running.')
        elif metadata['status'] == 'Failed':
            logger.info('This job has failed. Try to resubmit first.')
        elif metadata['status'] == 'Succeeded':
            outputs = metadata['outputs']
            # format outputs
            output_path = {}
            for key, value in outputs.items():
                output_name = key.split(".")[-1]
                if output_name == "unmapped_bams":
                    continue
                if isinstance(value, str):
                    output_path[output_name] = {value: os.path.join(sample_dst, os.path.basename(value))}
                elif isinstance(value, list):
                    this_group_output_dir = os.path.join(sample_dst, output_name)
                    if not os.path.exists(this_group_output_dir):
                        os.makedirs(this_group_output_dir)
                    output_path[output_name] = [{value[i]: os.path.join(this_group_output_dir, os.path.basename(value[i]))} for i in range(len(value))]
                    # for i in range(len(value)):
                    #     output_path[output_name+"_"+str(i+1)] = {value[i]: os.path.join(this_group_output_dir, os.path.basename(value[i]))}
                elif isinstance(value, float):
                    with open(os.path.join(sample_dst, output_name+".txt"), 'w') as file:
                        file.write(str(value))
            return output_path

    def delete_directories(self, metadata_list=None):
        if metadata_list is not None:
            for m in metadata_list:
                metadata = json.load(open(m, "r"))
                directory = metadata["workflowRoot"]
                os.rmdir(directory)

    def list(self, keyword=None, keyword_exclude=None):
        '''get a list of jobs on the server, their ids, names and statuses'''
        cmd = ['caper', 'list']

        if isinstance(self.hostname, str):
            cmd.extend(['--hostname', self.hostname])
        if isinstance(self.port, int):
            cmd.extend(['--port', str(self.port)])

        p = subprocess.run(cmd, capture_output=True)
        output, error = p.stdout.decode().strip(), p.stderr.decode().strip()
        jobs = output.split('\n')[1:]

        if keyword is not None:
            if isinstance(keyword, str):
                keyword = [keyword]
            for k in keyword:
                jobs = [j for j in jobs if k in j]
        if keyword_exclude is not None:
            if isinstance(keyword_exclude, str):
                keyword_exclude = [keyword_exclude]
            for k in keyword_exclude:
                jobs = [j for j in jobs if k not in j]

        job_ids = []
        statuses = []
        names = []
        for j in jobs:
            job_ids.append(j.split("\t")[0])
            statuses.append(j.split("\t")[1])
            names.append(j.split("\t")[3])
        
        return job_ids, statuses, names
    
    def count_list(self, keyword=None, keyword_exclude=None):
        job_ids, statuses, names = self.list(keyword=keyword,keyword_exclude=keyword_exclude)
        counts = Counter(statuses)
        logger.info('Out of {} workflows, {} succeeded, {} failed, {} still running'.format(len(statuses), counts['Succeeded'], counts['Failed'], counts['Running']))
    
    def debug(self, uid=None):
        '''Given workflow uid, get this workflow's error message and where the error happened, including the directory of the call.'''
        if uid is not None:
            cmd = ['caper', 'debug', uid]

            if isinstance(self.port, int):
                cmd.extend(['--port', str(self.port)])
            if isinstance(self.hostname, str):
                cmd.extend(['--hostname', self.hostname])

            p = subprocess.run(cmd, capture_output=True)
            output, error = p.stdout.decode().strip(), p.stderr.decode().strip()
            lines = [line.strip() for line in output.split('\n')]
            error_message = [line.replace('"message": ','') for line in lines if re.search('^"message"', line)]
            stderr = [line.split('=')[1] for line in lines if re.search('^STDERR=', line)]

            if len(error_message) != 0 and len(stderr) != 0:
                logger.info(' The error message is: {}\nCheck {} for details'.format(error_message[0], stderr[0]))
            # return error_message, stderr
    
    def _execute(self, cmd, silent=False):
        cmd = [str(c) for c in cmd]
        logger.info('Command: {}'.format(' '.join(cmd)))
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, bufsize=0)

        output = []
        for line in iter(process.stdout.readline, ''):
            if not silent:
                sys.stdout.write(line)
                sys.stdout.flush()
            output.append(line.strip())

        process.stdout.close()
        return(output)
    
    def _parse_run_output(self, output):
        line = [line for line in output if re.search(
            'Wrote metadata file.', line)][0]
        outfile = re.search('Wrote metadata file. (.*)', line).group(1)
        outdir = '/'.join(outfile.split('/')[0:-2])
        uid = outfile.split('/')[-2]
        return(outfile)

    def _parse_server_output(self, output):
        line = [line for line in output if re.search(
            "'id': \'(.*)\'\,", line)][0]
        uid = re.search("'id': \'(.*)\'\,", line).group(1)
        return(uid)

    def _process_job_outputs(self, outfile, log_path):
        if outfile[0:3] == 'gs:':
            self.metadata = json.load(hl.hadoop_open(outfile, 'r'))
        else:
            self.metadata = json.load(open(outfile, 'r'))
            
        json.dump(self.metadata, open(log_path, 'w', encoding='utf-8'), ensure_ascii = False, indent = 4)
        
        if os.path.exists(f'{self.filedir}/cromwell.out'):
            shutil.move(f'{self.filedir}/cromwell.out',
                os.path.splitext(log_path)[0] + '.out')

    def get_inputs(self, wdl):
        cmd = ["java", "-jar", self.womtool_jar, "inputs", wdl]
        p = subprocess.run(cmd, capture_output=True)
        output, error = p.stdout.decode().strip(), p.stderr.decode().strip()
        inputs = json.loads(output)
        #print(inputs)
        return(output)

    def get_outputs(self, metadata=None):
        metadata = self.f if metadata is None else metadata
        jdict = json.load(hl.hadoop_open(metadata, 'r'))
        return(jdict['outputs'])
    
    def get_metadata(self, uid=None):
        'server mode only'
        uid = uid if uid is not None else self.uid
        cmd = ['caper', 'metadata', uid]
        if isinstance(self.hostname, str):
            cmd.extend(['--hostname', self.hostname])
        if isinstance(self.port, int):
            cmd.extend(['--port', str(self.port)])
        p = subprocess.run(cmd, capture_output=True)
        output, error = p.stdout.decode().strip(), p.stderr.decode().strip()
        metadata = json.loads(output)
        return(metadata)
    
    def get_status(self, uid=None):
        metadata = self.get_metadata(uid)
        return(metadata['status'])

    def done(self, outfile, abspath=False):
        #----
        if self.metadata['status'] != 'Succeeded':
            raise Exception('Cromwell did not finish successfully.')

        self.fofn[self.current_run] = outfile
        self.current_run += 1
        if self._sig_popr:
            self._sig_popr = False
            self._setr.pop()
        
        #----
        with open(os.path.join(self.logdir, 'output.fofn'), 'w') as f:
            if abspath:
                fofn = [os.path.abspath(filename) for filename in self.fofn]
            else:
                fofn = self.fofn
            f.write('\n'.join(fofn) + '\n')

    def write_json(self, outfile, jdict):
        json.dump(jdict, open(outfile, 'w', encoding='utf-8'),
                  ensure_ascii=False, indent=4)

    def test_env(self):
        cmd = ['which', 'pip']
        p = subprocess.run(cmd, capture_output=True)
        _, error = p.stdout.decode().strip(), p.stderr.decode().strip()
        return(error)

    #---------------------------------------------------
    # hailrunner (properties)

    @property
    def r(self):
        return(self.current_run)
    
    @property
    def f(self):
        if len(self.fofn[self.current_run]) == 1:
            return(self.fofn[self.current_run][0])
        return(self.fofn[self.current_run])

    @property
    def popr(self):
        self._sig_popr = True
        return(self._setr[-1])

    @property
    def popf(self):
        self._sig_popr = True
        if len(self.fofn[self._setr[-1]]) == 1:
            return(self.fofn[self._setr[-1]][0])
        return(self.fofn[self._setr[-1]])

    # should be setter, but it sort of behaves differently
    def pushr(self):
        self._setr.append(self.current_run)

    @property
    def past(self):
        paths = [path.strip() for path in execute(['cat', os.path.join(self.logdir, 'output.fofn')],
                                                  return_out=True, silent=True)]
        return(list(enumerate(paths)))
    
    #---------------------------------------------------
    
    @property
    def ids(self):
        return([p.split('/')[-2] for p in self.fofn if p != ''])
    
    @property
    def outputs(self):
        return(self.metadata['outputs'])
    
    @property
    def outputs_list(self):
        return([v for k, v in self.metadata['outputs'].items()])
    
    @property
    def root(self):
        return(self.metadata['workflowRoot'])
