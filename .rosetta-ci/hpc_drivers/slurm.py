# -*- coding: utf-8 -*-
# :noTabs=true:

import os, sys, time, collections, math
import stat as stat_module


try:
    from .base import *

except ImportError: # workaround for B2 back-end's
    import imp
    imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) + '/base.py')  # A bit of Python magic here, what we trying to say is this: from base import *, but path to base is calculated from our source location  # from base import HPC_Driver, execute, NT


_T_slurm_array_job_template_ = '''\
#!/bin/bash
#
#SBATCH --job-name={name}
#SBATCH --output={log_dir}/.hpc.%x.%a.output
#
#SBATCH --time={time}:00
#SBATCH --mem-per-cpu={memory}M
#SBATCH --chdir={working_dir}
#
#SBATCH --array=1-{jobs_to_queue}

srun {executable} {arguments}
'''

_T_slurm_mpi_job_template_ = '''\
#!/bin/bash
#
#SBATCH --job-name={name}
#SBATCH --output={log_dir}/.hpc.%x.output
#
#SBATCH --time={time}:00
#SBATCH --mem-per-cpu={memory}M
#SBATCH --chdir={working_dir}
#
#SBATCH --ntasks={ntasks}

mpirun {executable} {arguments}
'''

class Slurm_HPC_Driver(HPC_Driver):
    def head_node_execute(self, message, command_line, *args, **kwargs):
        head_node = self.config['slurm'].get('head_node')

        command_line, host = (f"ssh {head_node} cd `pwd` '&& {command_line}'", head_node) if head_node else (command_line, 'localhost')
        return execute(f'Executiong on {host}: {message}' if message else '', command_line, *args, **kwargs)


    # NodeGroup = collections.namedtuple('NodeGroup', 'nodes cores')

    # @property
    # def mpi_topology(self):
    #     ''' return list of NodeGroup's
    #     '''
    #     pass


    # @property
    # def number_of_cpu_per_node(self): return int( self.config['condor']['mpi_cpu_per_node'] )

    # @property
    # def maximum_number_of_mpi_cpu(self):
    #     return self.number_of_cpu_per_node * int( self.config['condor']['mpi_maximum_number_of_nodes'] )


    # def complete(self, condor_job_id):
    #     ''' Return job completion status. Note that single hpc_job may contatin inner list of individual HPC jobs, True should be return if they all run in to completion.
    #     '''

    #     execute('Releasing condor jobs...', 'condor_release $USER', return_='tuple')

    #     s = execute('', 'condor_q $USER | grep $USER | grep {}'.format(condor_job_id), return_='output', terminate_on_failure=False).replace(' ', '').replace('\n', '')
    #     if s: return False

    #         # #setDaemonStatusAndPing('[Job #%s] Running... %s condor job(s) in queue...' % (self.id, len(s.split('\n') ) ) )
    #         # n_jobs = len(s.split('\n'))
    #         # s, o = execute('', 'condor_userprio -all | grep $USER@', return_='tuple')
    #         # if s == 0:
    #         #     jobs_running = o.split()
    #         #     jobs_running = 'XX' if len(jobs_running) < 4 else jobs_running[4]
    #         #     self.set_daemon_message("Waiting for condor to finish HPC jobs... [{} jobs in HPC-Queue, {} CPU's used]".format(n_jobs, jobs_running) )
    #         #     print "{} condor jobs in queue... Sleeping 32s...    \r".format(n_jobs),
    #         # sys.stdout.flush()
    #         # time.sleep(32)
    #     else:

    #         #self.tracer('Waiting for condor to finish the jobs... DONE')
    #         self.jobs.remove(condor_job_id)
    #         self.cpu_usage += self.get_condor_accumulated_usage()
    #         return True  # jobs already finished, we return empty list to prevent double counting of cpu_usage


    def complete(self, slurm_job_id):
        ''' Return True if job with given id is complete
        '''

        s = self.head_node_execute('', f'squeue -j {slurm_job_id} --noheader', return_='output', terminate_on_failure=False, silent=True)
        if s: return False
        else:
            #self.tracer('Waiting for condor to finish the jobs... DONE')
            self.jobs.remove(slurm_job_id)
            return True  # jobs already finished, we return empty list to prevent double counting of cpu_usage


    def cancel_job(self, slurm_job_id):
        self.head_node_execute(f'Slurm_HPC_Driver.canceling job {slurm_job_id}...', f'scancel {slurm_job_id}', terminate_on_failure=False)


    # def submit_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):
    #     print('submit_hpc_job is DEPRECATED and will be removed in near future, please use submit_serial_hpc_job  instead!')
    #     return self.submit_serial_hpc_job(name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory, time, block, shell_wrapper)


    def submit_serial_hpc_job(self, name, executable, arguments, working_dir, jobs_to_queue, log_dir, memory=512, time=12, block=True, shell_wrapper=False):

        arguments = arguments.format(process='%a') # %a is SLURM array index
        time = int( math.ceil(time*60) )

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + f'/hpc.{name}.shell_wrapper.sh')
            with open(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        slurm_file = working_dir + f'/.hpc.{name}.slurm'

        with open(slurm_file, 'w') as f: f.write( _T_slurm_array_job_template_.format( **vars() ) )


        slurm_job_id = self.head_node_execute('Submitting SLURM array job...', f'cd {self.working_dir} && sbatch {slurm_file}',
                                              tracer=self.tracer, return_='output'
                                              ).split()[-1] # expecting something like `Submitted batch job 6122` in output


        self.jobs.append(slurm_job_id)

        if block:
            self.wait_until_complete( [slurm_job_id] )
            return None

        else: return slurm_job_id





    def submit_mpi_hpc_job(self, name, executable, arguments, working_dir, log_dir, ntasks, memory=512, time=12, block=True, shell_wrapper=False):
        ''' submit jobs as MPI job
        '''
        arguments = arguments.format(process='0')
        time = int( math.ceil(time*60) )

        if shell_wrapper:
            shell_wrapper_sh = os.path.abspath(self.working_dir + f'/hpc.{name}.shell_wrapper.sh')
            with open(shell_wrapper_sh, 'w') as f: f.write('#!/bin/bash\n{} {}\n'.format(executable, arguments));  os.fchmod(f.fileno(), stat.S_IEXEC | stat.S_IREAD | stat.S_IWRITE)
            executable, arguments = shell_wrapper_sh, ''

        slurm_file = working_dir + f'/.hpc.{name}.slurm'

        with open(slurm_file, 'w') as f: f.write( _T_slurm_mpi_job_template_.format( **vars() ) )

        slurm_job_id = self.head_node_execute('Submitting SLURM mpi job...', f'cd {self.working_dir} && sbatch {slurm_file}',
                                              tracer=self.tracer, return_='output'
                                              ).split()[-1] # expecting something like `Submitted batch job 6122` in output

        self.jobs.append(slurm_job_id)

        if block:
            self.wait_until_complete( [slurm_job_id] )
            return None

        else: return slurm_job_id
