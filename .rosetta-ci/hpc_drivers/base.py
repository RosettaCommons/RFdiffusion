# -*- coding: utf-8 -*-
# :noTabs=true:

import os, sys, subprocess, stat
import time as time_module
import signal as signal_module

class NT:  # named tuple
    def __init__(self, **entries): self.__dict__.update(entries)
    def __repr__(self):
        r = 'NT: |'
        for i in dir(self):
            if not i.startswith('__') and not isinstance(getattr(self, i), types.MethodType): r += '{} --> {}, '.format(i, getattr(self, i))
        return r[:-2]+'|'



class HPC_Exception(Exception):
    def __init__(self, value): self.value = value
    def __str__(self): return self.value



def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False, tracer=print):
    if not silent: tracer(message);  tracer(command_line); sys.stdout.flush();
    while True:

        p = subprocess.Popen(command_line, bufsize=0, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, errors = p.communicate()

        output = output + errors

        output = output.decode(encoding="utf-8", errors="replace")

        exit_code = p.returncode

        if exit_code  and  not (silent or silence_output): tracer(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        tracer( "Error while executing {}: {}\n".format(message, output) )
        tracer("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if return_ == 'tuple': return(exit_code, output)

    if exit_code and terminate_on_failure:
        tracer("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else: print("\nEncounter error while executing: " + command_line + '\n' + output); sys.exit(1)

    if return_ == 'output': return output
    else: return False


def Sleep(time_, message, dict_={}):
    ''' Fancy sleep function '''
    len_ = 0
    for i in range(time_, 0, -1):
        #print "Waiting for a new revision:%s... Sleeping...%d     \r" % (sc.revision, i),
        msg = message.format( **dict(dict_, time_left=i) )
        print( msg, end='' )
        len_ = max(len_, len(msg))
        sys.stdout.flush()
        time_module.sleep(1)

    print( ' '*len_ + '\r',  end='' ) # erazing sleep message


# Abstract class for HPC job submission
class HPC_Driver:
    def __init__(self, working_dir, config, tracer=lambda x:None, set_daemon_message=lambda x:None):
        self.working_dir = working_dir
        self.config = config
        self.cpu_usage  = 0.0  # cummulative cpu usage in hours
        self.tracer     = tracer
        self.set_daemon_message = set_daemon_message

        self.cpu_count = self.config['cpu_count'] if type(config) == dict else self.config.getint('DEFAULT', 'cpu_count')

        self.jobs = []  # list of all jobs currently running by this driver, Job class is driver depended, could be just int or something more complex

        self.install_signal_handler()


    def __del__(self):
        self.remove_signal_handler()


    def execute(self, executable, arguments, working_dir, log_dir=None, name='_no_name_', memory=256, time=24, shell_wrapper=False, block=True):
        ''' Execute given command line on HPC cluster, must accumulate cpu hours in self.cpu_usage '''
        if log_dir==None: log_dir=self.working_dir

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + '/hpc.{}.shell_wrapper.sh'.format(name))
            with file(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        return self.submit_serial_hpc_job(name=name, executable=executable, arguments=arguments, working_dir=working_dir, log_dir=log_dir, jobs_to_queue=1, memory=memory, time=time, block=block, shell_wrapper=shell_wrapper)



    @property
    def number_of_cpu_per_node(self):
        must_be_implemented_in_inherited_classes

    @property
    def maximum_number_of_mpi_cpu(self):
        must_be_implemented_in_inherited_classes


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        print('submit_hpc_job is DEPRECATED and will be removed in near future, please use submit_serial_hpc_job  instead!')
        must_be_implemented_in_inherited_classes


    def submit_serial_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        must_be_implemented_in_inherited_classes


    def submit_mpi_hpc_job(self, name, executable, arguments, working_dir, log_dir, memory=512, time=12, block=True, process_coefficient="1", requested_nodes=1, requested_processes_per_node=1):
        ''' submit jobs as MPI job
            process_coefficient should be string representing fraction of process to launch on each node, for example '3 / 4' will start only 75% of MPI process's on each node
        '''
        must_be_implemented_in_inherited_classes


    def cancel_all_jobs(self):
        ''' Cancel all HPC jobs known to this driver, use this as signal handler for script termination '''
        for j in self.jobs: self.cancel_job(j)

    def block_until(self, silent, fn, *args, **kwargs):
        '''
        **fn must have the driver as the first argument**
        example:
        def fn(driver):
            jobs = list(driver.jobs)
            jobs = [job for job in jobs if not driver.complete(job)]
            if len(jobs) <= 8:
                return False # stops sleeping
            return True # continues sleeping

        for x in range(100):
            hpc_driver.submit_hpc_job(...)
            hpc_driver.block_until(False, fn)
        '''
        while fn(self, *args, **kwargs):
            sys.stdout.flush()
            time_module.sleep(60)
            if not silent:
                Sleep(1, '"Waiting for HPC job(s) to finish, sleeping {time_left}s\r')

    def wait_until_complete(self, jobs=None, callback=None, silent=False):
        ''' Helper function, wait until given jobs list is finished, if no argument is given waits until all jobs known by driver is finished '''
        jobs = jobs if jobs else self.jobs

        while jobs:
            for j in jobs[:]:
                if self.complete(j): jobs.remove(j)

            if jobs:
                #total_cpu_queued  = sum( [j.jobs_queued  for j in jobs] )
                #total_cpu_running = sum( [j.cpu_running for j in jobs] )
                #self.set_daemon_message("Waiting for HPC job(s) to finish... [{} process(es) in queue, {} process(es) running]".format(total_cpu_queued, total_cpu_running) )
                #self.tracer("Waiting for HPC job(s) [{} process(es) in queue, {} process(es) running]...  \r".format(total_cpu_queued, total_cpu_running), end='')
                #print "Waiting for {} HPC jobs to finish... [{} jobs in queue, {} jobs running]... Sleeping 32s...     \r".format(total_cpu_queued, cpu_queued+cpu_running, cpu_running),

                self.set_daemon_message("Waiting for HPC {} job(s) to finish...".format( len(jobs) ) )
                #self.tracer("Waiting for HPC {} job(s) to finish...".format( len(jobs) ) )

                sys.stdout.flush()

                if callback: callback()

                if silent: time_module.sleep(64*1)
                else: Sleep(64, '"Waiting for HPC {n_jobs} job(s) to finish, sleeping {time_left}s    \r', dict(n_jobs=len(jobs)))



    _signals_ = [signal_module.SIGINT, signal_module.SIGTERM, signal_module.SIGABRT]
    def install_signal_handler(self):
        def signal_handler(signal_, frame):
            self.tracer('Recieved signal:{}... Canceling HPC jobs...'.format(signal_) )
            self.cancel_all_jobs()
            self.set_daemon_message( 'Remote daemon got terminated with signal:{}'.format(signal_) )
            sys.exit(1)

        for s in self._signals_: signal_module.signal(s, signal_handler)


    def remove_signal_handler(self):  # do we really need this???
        try:
            for s in self._signals_: signal_module.signal(s, signal_module.SIG_DFL)
            #print('remove_signal_handler: done!')

        except TypeError:
            #print('remove_signal_handler: interpreted terminating, skipping remove_signal_handler...')
            pass


    def cancel_job(self, job_id):
        must_be_implemented_in_inherited_classes


    def complete(self, job_id):
        ''' Return job completion status. Return True if job complered and False otherwise
        '''
        must_be_implemented_in_inherited_classes
