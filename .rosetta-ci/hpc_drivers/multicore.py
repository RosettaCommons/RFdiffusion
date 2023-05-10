# -*- coding: utf-8 -*-
# :noTabs=true:

import time as time_module
import codecs
import signal

import os, sys

try:
    from .base import *

except ImportError: # workaround for B2 back-end's
    import imp
    imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


class MultiCore_HPC_Driver(HPC_Driver):

    class JobID:
        def __init__(self, pids=None):
            self.pids = pids if pids else []


        def __bool__(self): return bool(self.pids)


        def __len__(self): return len(self.pids)


        def add_pid(self, pid): self.pids.append(pid)


        def remove_completed_pids(self):
            for pid in self.pids[:]:
                try:
                    r = os.waitpid(pid, os.WNOHANG)
                    if r == (pid, 0): self.pids.remove(pid)  # process have ended without error
                    elif r[0] == pid :  # process ended but with error, special case we will have to wait for all process to terminate and call system exit.
                        #self.cancel_job()
                        #sys.exit(1)
                        self.pids.remove(pid)
                        print('ERROR: Some of the HPC jobs terminated abnormally! Please see HPC logs for details.')

                except ChildProcessError: self.pids.remove(pid)


        def cancel(self):
            for pid in self.pids:
                try:
                    os.killpg(os.getpgid(pid), signal.SIGKILL)
                except ChildProcessError: pass

            self.pids = []



    def __init__(self, *args, **kwds):
        HPC_Driver.__init__(self, *args, **kwds)
        #print(f'MultiCore_HPC_Driver: cpu_count: {self.cpu_count}')


    def remove_completed_jobs(self):
        for job in self.jobs[:]:  # Need to make a copy so we don't modify a list we're iterating over
            job.remove_completed_pids()
            if not job: self.jobs.remove(job)


    @property
    def process_count(self):
        ''' return number of processes that currently ran by this driver instance
        '''
        return sum( map(len, self.jobs) )


    def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        print('submit_hpc_job is DEPRECATED and will be removed in near future, please use submit_serial_hpc_job  instead!')
        return self.submit_serial_hpc_job(name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory, time, block, shell_wrapper)


    def submit_serial_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
        cpu_usage = -time_module.time()/60./60.

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + f'/hpc.{name}.shell_wrapper.sh')
            with open(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        def mfork():
            ''' Check if number of child process is below cpu_count. And if it is - fork the new pocees and return its pid.
            '''
            while self.process_count >= self.cpu_count:
                self.remove_completed_jobs()
                if self.process_count >= self.cpu_count: time_module.sleep(.5)

            sys.stdout.flush()
            pid = os.fork()
            # appending at caller level insted  if pid: self.jobs.append(pid) # We are parent!
            return pid

        current_job = self.JobID()
        process = 0
        for i in range(jobs_to_queue):

            pid = mfork()
            if not pid: # we are child process
                command_line = 'cd {} && {} {}'.format(working_dir, executable, arguments.format(process=process) )
                exit_code, log = execute('Running job {}.{}...'.format(name, i), command_line, tracer=self.tracer, return_='tuple')
                with codecs.open(log_dir+'/.hpc.{name}.{i:02d}.log'.format(**vars()), 'w', encoding='utf-8', errors='replace') as f:
                    f.write(command_line+'\n'+log)
                    if exit_code:
                        error_report = f'\n\n{command_line}\nERROR: PROCESS {name}.{i:02d} TERMINATED WITH NON-ZERO-EXIT-CODE {exit_code}!\n'
                        f.write(error_report)
                        print(log, error_report)

                sys.exit(0)

            else: # we are parent!
                current_job.add_pid(pid)
                # Need to potentially re-add to list, as remove_completed_jobs() might trim it.
                if current_job not in self.jobs: self.jobs.append(current_job)

            process += 1

        if block:
            #for p in all_queued_jobs: os.waitpid(p, 0)  # waiting for all child process to termintate...

            self.wait_until_complete(current_job)
            self.remove_completed_jobs()

            cpu_usage += time_module.time()/60./60.
            self.cpu_usage += cpu_usage * jobs_to_queue  # approximation...

            current_job = self.JobID()

        return current_job


    @property
    def number_of_cpu_per_node(self): return self.cpu_count


    @property
    def maximum_number_of_mpi_cpu(self): return self.cpu_count


    def submit_mpi_hpc_job(self, name, executable, arguments, working_dir, log_dir, memory=512, time=12, block=True, process_coefficient="1", requested_nodes=1, requested_processes_per_node=1):

        if requested_nodes > 1:
            print( "WARNING: " + str( requested_nodes ) + " nodes were requested, but we're running locally, so only 1 node will be used." )

        if requested_processes_per_node > self.cpu_count:
            print( "WARNING: " + str(requested_processes_per_node) + " processes were requested, but I only have " + str(self.cpu_count) + " CPUs.  Will launch " + str(self.cpu_count) + " processes."  )
        actual_processes = min( requested_processes_per_node, self.cpu_count )

        cpu_usage = -time_module.time()/60./60.

        arguments = arguments.format(process=0)

        command_line = f'cd {working_dir} && mpirun -np {actual_processes} {executable} {arguments}'
        log = execute(f'Running job {name}...', command_line, tracer=self.tracer, return_='output')
        with codecs.open(log_dir+'/.hpc.{name}.log'.format(**vars()), 'w', encoding='utf-8', errors='replace') as f: f.write(command_line+'\n'+log)

        cpu_usage += time_module.time()/60./60.
        self.cpu_usage += cpu_usage * actual_processes  # approximation...

        # return None - we do not return anything from this version of submit which imply returning None which in turn will be treated as job-id for already finished job


    def complete(self, job_id):
        ''' Return job completion status. Return True if job completed and False otherwise
        '''
        self.remove_completed_jobs()
        return job_id not in self.jobs


    def cancel_job(self, job):
        job.cancel();
        if job in self.jobs:
            self.jobs.remove(job)


    def __repr__(self):
        return 'MultiCore_HPC_Driver<>'
