import unittest
import subprocess
import glob
import datetime
import os
import torch
from shutil import copyfile
from rfdiffusion.inference import utils as iu
from rfdiffusion.util import calc_rmsd
import sys

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
        for bash_file in glob.glob(f"{self.out_f}/*.sh"):
            print(f"Running {os.path.basename(bash_file)}")
            subprocess.run(["bash", bash_file], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    def test_commands(self):
        """
        Runs all the commands in the test_f folder
        """
        reference=f'{script_dir}/reference_outputs'
        os.makedirs(reference, exist_ok=True)
        test_files=glob.glob(f"{self.out_f}/example_outputs/*pdb")
        # first check that we have the right number of outputs
        self.assertEqual(len(test_files), len(glob.glob(f"{self.out_f}/*.sh"))), "One or more of the example commands didn't produce an output (check the example command is formatted correctly)"
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
                    except AssertionError as e:
                        result.addFailure(self, e)
                        print(f"Subtest {test_file} failed")
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

if __name__ == "__main__":
    unittest.main()
