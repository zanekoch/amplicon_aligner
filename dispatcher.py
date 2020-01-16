import logging
import aligner

def _submit_job(job_submitter, command, job_parms, waitfor_id=None, hold=False, notify=False):
    import subprocess
    import re
    import os
    # TODO(jtravis): remove unused output variable
    output = jobid = None
    logging.info("command = %s" % command)
    args = job_parms["args"]
    if job_submitter == "PBS":
        waitfor = ""
        if waitfor_id and waitfor_id[0]:
            dependency_string = waitfor_id[1] if len(waitfor_id) > 1 else 'afterok'
            waitfor = "-W depend=%s:%s" % (dependency_string, waitfor_id[0])
        queue = ""
        if job_parms["queue"]:
            queue = "-q %s" % job_parms["queue"]
        if hold:
            args += " -h"
        if notify:
            args += " -m e"
        submit_command = "qsub -V -d \'%s\' -w \'%s\' -l ncpus=%s,mem=%sgb,walltime=%s:00:00 -m a -N \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'],
            job_parms['walltime'], job_parms['name'], waitfor, queue, args)
        logging.debug("submit_command = %s", submit_command)
        output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
        logging.debug("output = %s" % output)
        job_match = re.search('^(\d+)\..*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SLURM":
        waitfor = ""
        if waitfor_id:
            dependency_string = waitfor_id[1] if len(waitfor_id) > 1 else 'afterok'
            waitfor = "-d %s:%s" % (dependency_string, waitfor_id[0])
        queue = ""
        if job_parms["queue"]:
            queue = "-p %s" % job_parms["queue"]
        if hold:
            args += " -H"
        if notify:
            args += " --mail-type=END"
        #submit_command = "sbatch -D \'%s\' -c%s --mem=%s000 --time=%s:00:00 --mail-type=FAIL -J \'%s\' %s %s %s" % (
        #    job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'], job_parms['walltime'],
        #    job_parms['name'], waitfor, queue, args)
        submit_command = "sbatch -D \'%s\' -c%s --mem=%s000 --mail-type=FAIL -J \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms['num_cpus'], job_parms['mem_requested'],
            job_parms['name'], waitfor, queue, args)
        logging.debug("submit_command = %s" % submit_command)
        output = subprocess.getoutput("%s --wrap=\"%s\"" % (submit_command, command))
        logging.debug("output = %s" % output)
        job_match = re.search('^Submitted batch job (\d+)$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    elif job_submitter == "SGE":
        waitfor = ""
        if waitfor_id:
            waitfor = "-hold_jid %s" % (re.sub(":", ",", waitfor_id[0]))
        queue = ""
        if job_parms["queue"]:
            queue = "-q %s" % job_parms["queue"]
        if hold:
            args += " -h"
        if notify:
            args += " -m e"
        mem_needed = float(job_parms['mem_requested']) * 1024 * 1024
        # Apparently the number of processors a job uses is controlled by the queue it is running on in SGE, so there is no way to request a specific number of CPUs??
        submit_command = "qsub -V -cwd \'%s\' -wd \'%s\' -l h_data=%sgb,h_rt=%s:00:00 -m a -N \'%s\' %s %s %s" % (
            job_parms["work_dir"], job_parms["work_dir"], mem_needed, job_parms['walltime'], job_parms['name'], waitfor,
            queue, args)
        logging.debug("submit_command = %s", submit_command)
        output = subprocess.getoutput("echo \"%s\" | %s - " % (command, submit_command))
        logging.debug("output = %s" % output)
        job_match = re.search('^(\d+)\..*$', output)
        if job_match:
            jobid = job_match.group(1)
        else:
            logging.warning("Job not submitted!!")
            print("WARNING: Job not submitted: %s" % output)
    else:
        command = re.sub('\n', '; ', command)
        work_dir = job_parms['work_dir']
        dependency_check = ""
        if waitfor_id:
            if re.search(":", waitfor_id[0]):
                pid_filename = os.path.join(work_dir, "%s_dependent_pids" % job_parms['name'])
                pid_file = open(pid_filename, 'w')
                pids = waitfor_id[0].split(":")
                pid_file.write("\n".join(pids))
                pid_file.close()
                dependency_check = "while [ -s %s ]; do sleep 600; for pid in `cat %s`; do kill -0 \"$pid\" 2>/dev/null || sed -i \"/^$pid$/d\" %s; done; done; rm %s; " % (
                    pid_filename, pid_filename, pid_filename, pid_filename)
            else:
                pid = waitfor_id[0]
                dependency_check = "while kill -0 %s; do sleep 300; done; " % pid
        # cpu_check = "grep 'cpu ' /proc/stat | awk '{usage=($2+$4)/($2+$4+$5)} END {print 1-usage}'"
        total_mem = subprocess.getoutput("free -g | grep Mem: | awk '{ print $2 }'")
        mem_requested = job_parms['mem_requested']
        if mem_requested > total_mem:
            mem_requested = total_mem
        print("total_mem = %s" % total_mem)
        mem_needed = float(mem_requested) * 750
        mem_check = "while [ `free -m | grep cache: | awk '{ print $4 }'` -lt %d ]; do sleep 300; done; " % mem_needed
        submit_command = "%s%s%s" % (dependency_check, mem_check, command)
        logging.debug("submit_command = %s" % submit_command)
        output_log = os.open(os.path.join(work_dir, "%s.out" % job_parms['name']), os.O_WRONLY | os.O_CREAT)
        proc = subprocess.Popen(submit_command, stderr=subprocess.STDOUT, stdout=output_log, shell=True, cwd=work_dir)
        jobid = str(proc.pid)
    logging.info("jobid = %s" % jobid)
    return jobid


def startAlignment(read_tuple, ref_filename, out_dir, job_manager, hash_length):
    #create a job for each read1-read2 file pair that puts these reads into memory
    job_params = {'queue':'', 'mem_requested':8, 'num_cpus':2, 'walltime':24, 'args':''}
    #read_tuple[0] is sample name
    job_params['name'] = "retrieveReads_%s" % read_tuple[0]
    job_params['work_dir'] = out_dir
    read1 = read_tuple[1][0]
    if len(read_tuple[1]) > 1:
        read2 = read_tuple[1][1]
    else:
        read2 = ''
    command = "python /scratch/zkoch/amplicon_aligner/aligner.py -s %s -r1 %s -r2 %s -f %s -o %s -j %s --hash-length %d" % (read_tuple[0], read1, read2, ref_filename, out_dir, job_manager, hash_length)
    jobid = _submit_job(job_manager, command, job_params)

    return jobid
