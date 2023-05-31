import unittest
import subprocess
import glob
import datetime
import os
import torch
from shutil import copyfile
from rfdiffusion.inference import utils as iu
from rfdiffusion.util import calc_rmsd
import sys, json

script_dir = os.path.dirname(os.path.abspath(__file__))

class TestSubmissionCommands(unittest.TestCase):
    """
    Test harness for checking that commands in the examples folder,
    when run in deterministic mode, produce the same output as the
    reference outputs.
    Requirements:
        - example command must be written on a single line
        - outputs must be written to example_outputs folder
        - needs to be run on the same hardware as the reference outputs (A100 GPU)
    For speed, we only run the first 2 steps of diffusion, and set inference.num_designs=1
    This means that outputs DO NOT look like proteins, but we can still check that the
    outputs are the same as the reference outputs.
    """

    def setUp(self):
        """
        Grabs files from the examples folder
        """
        submissions = glob.glob(f"{script_dir}/../examples/*.sh")
        # get datetime for output folder, in YYYY_MM_DD_HH_MM_SS format
        now = datetime.datetime.now()
        now = now.strftime("%Y_%m_%d_%H_%M_%S")
        self.out_f = f"{script_dir}/tests_{now}"
        os.mkdir(self.out_f)

        # Make sure we have access to all the relevant files
        exclude_dirs = ["outputs", "example_outputs"]
        for filename in os.listdir(f"{script_dir}/../examples"):
            if filename not in exclude_dirs and not os.path.islink(os.path.join(script_dir, filename)) and os.path.isdir(os.path.join(f'{script_dir}/../examples', filename)):
                os.symlink(os.path.join(f'{script_dir}/../examples', filename), os.path.join(script_dir, filename))

        for submission in submissions:
            self._write_command(submission, self.out_f)

        print(f"Running commands in {self.out_f}, two steps of diffusion, deterministic=True")

        self.results = {}

        for bash_file in sorted( glob.glob(f"{self.out_f}/*.sh"), reverse=False):
            test_name = os.path.basename(bash_file)[:-len('.sh')]
            res, output = execute(f"Running {test_name}", f'bash {bash_file}', return_='tuple', add_message_and_command_line_to_output=True)

            self.results[test_name] = dict(
                state = 'failed' if res else 'passed',
                log = output,
            )

            #subprocess.run(["bash", bash_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            #subprocess.run(["bash", bash_file])


    def test_commands(self):
        """
        Runs all the commands in the test_f folder
        """
        reference=f'{script_dir}/reference_outputs'
        os.makedirs(reference, exist_ok=True)
        test_files=glob.glob(f"{self.out_f}/example_outputs/*pdb")
        print(f'{self.out_f=} {test_files=}')

        # first check that we have the right number of outputs
        #self.assertEqual(len(test_files), len(glob.glob(f"{self.out_f}/*.sh"))), "One or more of the example commands didn't produce an output (check the example command is formatted correctly)"

        result = self.defaultTestResult()
        for test_file in test_files:
            with self.subTest(test_file=test_file):
                test_pdb=iu.parse_pdb(test_file)
                if not os.path.exists(f"{reference}/{os.path.basename(test_file)}"):
                    copyfile(test_file, f"{reference}/{os.path.basename(test_file)}")
                    print(f"Created reference file {reference}/{os.path.basename(test_file)}")
                else:
                    ref_pdb=iu.parse_pdb(f"{reference}/{os.path.basename(test_file)}")
                    rmsd=calc_rmsd(test_pdb['xyz'][:,:3].reshape(-1,3), ref_pdb['xyz'][:,:3].reshape(-1,3))[0]
                    try:
                        self.assertAlmostEqual(rmsd, 0, 2)
                        result.addSuccess(self)
                        print(f"Subtest {test_file} passed")

                        state = 'passed'
                        log = f'Subtest {test_file} passed'

                    except AssertionError as e:
                        result.addFailure(self, e)
                        print(f"Subtest {test_file} failed")

                        state = 'failed'
                        log = f'Subtest {test_file} failed:\n{e!r}'

                    self.results[ 'pdb-diff.' + test_file.rpartition('/')[-1] ] = dict(state = state, log = log)

        with open('.results.json', 'w') as f: json.dump(self.results, f, sort_keys=True, indent=2)

        self.assertTrue(result.wasSuccessful(), "One or more subtests failed")


    def _write_command(self, bash_file, test_f) -> None:
        """
        Takes a bash file from the examples folder, and writes
        a version of it to the test_f folder.
        It appends to the python command the following arguments:
            inference.deterministic=True
            if partial_T is in the command, it grabs partial T and sets:
                inference.final_step=partial_T-2
            else:
                inference.final_step=48
        """
        out_lines=[]
        with open(bash_file, "r") as f:
            lines = f.readlines()
            for line in lines:
                if not (line.startswith("python") or line.startswith("../")):
                    out_lines.append(line)
                else:
                    command = line.strip()
            if not command.startswith("python"):
                command = f'python {command}'
        # get the partial_T
        if "partial_T" in command:
            final_step = int(command.split("partial_T=")[1].split(" ")[0]) - 2
        else:
            final_step = 48

        output_command = f"{command} inference.deterministic=True inference.final_step={final_step}"
        # replace inference.num_designs with 1
        if "inference.num_designs=" in output_command:
            output_command = f'{output_command.split("inference.num_designs=")[0]}inference.num_designs=1 {" ".join(output_command.split("inference.num_designs=")[1].split(" ")[1:])}'
        else:
            output_command = f'{output_command} inference.num_designs=1'
        # replace 'example_outputs' with f'{self.out_f}/example_outputs'
        output_command = f'{output_command.split("example_outputs")[0]}{self.out_f}/example_outputs{output_command.split("example_outputs")[1]}'


        # write the new command
        with open(f"{test_f}/{os.path.basename(bash_file)}", "w") as f:
            for line in out_lines:
                f.write(line)
            f.write(output_command)



def execute_through_pty(command_line):
    import pty, select

    if sys.platform == "darwin":

        master, slave = pty.openpty()
        p = subprocess.Popen(command_line, shell=True, stdout=slave, stdin=slave,
                             stderr=subprocess.STDOUT, close_fds=True)

        buffer = []
        while True:
            try:
                if select.select([master], [], [], 0.2)[0]:  # has something to read
                    data = os.read(master, 1 << 22)
                    if data: buffer.append(data)

                elif (p.poll() is not None)  and  (not select.select([master], [], [], 0.2)[0] ): break  # process is finished and output buffer if fully read

            except OSError: break  # OSError will be raised when child process close PTY descriptior

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')

        os.close(master)
        os.close(slave)

        p.wait()
        exit_code = p.returncode

        '''
        buffer = []
        while True:
            if select.select([master], [], [], 0.2)[0]:  # has something to read
                data = os.read(master, 1 << 22)
                if data: buffer.append(data)
                # else: break  # # EOF - well, technically process _should_ be finished here...

            # elif time.sleep(1) or (p.poll() is not None): # process is finished (sleep here is intentional to trigger race condition, see solution for this on the next few lines)
            #     assert not select.select([master], [], [], 0.2)[0]  # should be nothing left to read...
            #     break

            elif (p.poll() is not None)  and  (not select.select([master], [], [], 0.2)[0] ): break  # process is finished and output buffer if fully read

        assert not select.select([master], [], [], 0.2)[0]  # should be nothing left to read...

        os.close(slave)
        os.close(master)

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')
        exit_code = p.returncode
        '''

    else:

        master, slave = pty.openpty()
        p = subprocess.Popen(command_line, shell=True, stdout=slave, stdin=slave,
                             stderr=subprocess.STDOUT, close_fds=True)

        os.close(slave)

        buffer = []
        while True:
            try:
                data = os.read(master, 1 << 22)
                if data: buffer.append(data)
            except OSError: break  # OSError will be raised when child process close PTY descriptior

        output = b''.join(buffer).decode(encoding='utf-8', errors='backslashreplace')

        os.close(master)

        p.wait()
        exit_code = p.returncode

    return exit_code, output



def execute(message, command_line, return_='status', until_successes=False, terminate_on_failure=True, silent=False, silence_output=False, silence_output_on_errors=False, add_message_and_command_line_to_output=False):
    if not silent: print(message);  print(command_line); sys.stdout.flush();
    while True:

        #exit_code, output = execute_through_subprocess(command_line)
        #exit_code, output = execute_through_pexpect(command_line)
        exit_code, output = execute_through_pty(command_line)

        if (exit_code  and  not silence_output_on_errors) or  not (silent or silence_output): print(output); sys.stdout.flush();

        if exit_code and until_successes: pass  # Thats right - redability COUNT!
        else: break

        print( "Error while executing {}: {}\n".format(message, output) )
        print("Sleeping 60s... then I will retry...")
        sys.stdout.flush();
        time.sleep(60)

    if add_message_and_command_line_to_output: output = message + '\nCommand line: ' + command_line + '\n' + output

    if return_ == 'tuple'  or  return_ == tuple: return(exit_code, output)

    if exit_code and terminate_on_failure:
        print("\nEncounter error while executing: " + command_line)
        if return_==True: return True
        else:
            print('\nEncounter error while executing: ' + command_line + '\n' + output);
            raise BenchmarkError('\nEncounter error while executing: ' + command_line + '\n' + output)

    if return_ == 'output': return output
    else: return exit_code


if __name__ == "__main__":
    unittest.main()
