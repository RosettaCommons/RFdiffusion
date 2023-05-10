#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rfd.py
## @brief  main test files for RFdiffusion
## @author Sergey Lyskov


import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'


def run_main_test_suite(repository_root, working_dir, platform, config):
    python_environment = local_python_install(platform, config)
    #ve = setup_persistent_python_virtual_environment(python_environment, packages='numpy torch omegaconf scipy opt_einsum dgl')
    ve = setup_python_virtual_environment(working_dir+'/.ve', python_environment, packages='numpy torch omegaconf scipy opt_einsum dgl e3nn icecream pyrsistent wandb pynvml decorator jedi hydra-core')
    execute('Installing local se3-transformer package...', f'cd {repository_root}/env/SE3Transformer && {ve.bin}/pip3 install -e .')

    print(f'{ve.python=}')

    #with tempfile.TemporaryDirectory() as tmpdirname:
    #    print('created temporary directory', tmpdirname)


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if test == '': return run_main_test_suite(repository_root, working_dir, platform, config)
    else: raise BenchmarkError('Unknow scripts test: {}!'.format(test))
