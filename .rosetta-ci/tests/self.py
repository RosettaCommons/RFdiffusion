#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   dummy.py
## @brief  self-test and debug-aids tests
## @author Sergey Lyskov

import os, os.path, shutil, re, string
import json

import random

import imp
imp.load_source(__name__, '/'.join(__file__.split('/')[:-1]) +  '/__init__.py')  # A bit of Python magic here, what we trying to say is this: from __init__ import *, but init is calculated from file location

_api_version_ = '1.0'


def run_state_test(repository_root, working_dir, platform, config):
    revision_id = config['revision']
    states = (_S_passed_, _S_failed_, _S_build_failed_, _S_script_failed_)
    state = states[revision_id % len(states)]

    return {_StateKey_ : state,  _ResultsKey_ : {},  _LogKey_ : f'run_state_test: setting test state to {state!r}...' }


sub_test_description_template = '''\
# subtests_test test suite
These sub-test description is generated for 3/4 of sub-tests

-----
### {name}
The warm time, had already disappeared like dust. Broken rain, fragment of light shadow, bring more pain to my heart...
-----
'''

def run_subtests_test(repository_root, working_dir, platform, config):
    tests = {}
    for i in range(16):
        name = f's-{i:02}'
        log = ('x'*63 + '\n') * 16 * 256 * i
        s = i % 3
        if   s == 0: state = _S_passed_
        elif s == 1: state = _S_failed_
        else:        state = _S_script_failed_

        if i % 4:
            os.mkdir( f'{working_dir}/{name}' )
            with open(f'{working_dir}/{name}/description.md', 'w') as f: f.write( sub_test_description_template.format(**vars()) )

            with open( f'{working_dir}/{name}/fantome.txt', 'w') as f: f.write('No one wants to hear the sequel to a fairytale\n')

        tests[name] = { _StateKey_ : state, _LogKey_   : log, }

    test_log = ('*'*63 + '\n') * 16 * 1024 * 16
    return {_StateKey_ : _S_failed_,  _ResultsKey_ : {_TestsKey_: tests},  _LogKey_ : test_log }


def run_regression_test(repository_root, working_dir, platform, config):
    const     = 'const'
    volatile  = 'volatile'
    new       = ''.join( random.sample( string.ascii_letters + string.digits, 8) )
    oversized = 'oversized'

    sub_tests = [const, volatile, new]

    const_dir = working_dir + '/' + const
    os.mkdir(const_dir)
    with open(const_dir + '/const_data', 'w') as f: f.write( '\n'.join( (str(i) for i in range(32) ) ) )

    volatile_dir = working_dir + '/' + volatile
    os.mkdir(volatile_dir)
    with open(volatile_dir + '/const_data', 'w') as f: f.write( '\n'.join( (str(i) for i in range(32, 64) ) ) )
    with open(volatile_dir + '/volatile_data', 'w') as f: f.write( '\n'.join( ( ''.join(random.sample( string.ascii_letters + string.digits, 8) ) for i in range(32) ) ) )

    new_dir = working_dir + '/' + new
    os.mkdir(new_dir)
    with open(new_dir + '/data', 'w') as f: f.write( '\n'.join( (str(i) for i in range(64)) ) )


    new_dir = working_dir + '/' + oversized
    os.mkdir(new_dir)
    with open(new_dir + '/large', 'w') as f: f.write( ('x'*63 + '\n')*16*1024*256 +'extra')

    return {_StateKey_ : _S_queued_for_comparison_,  _ResultsKey_ : {},  _LogKey_ : f'sub-tests: {sub_tests!r}' }



def run_release_test(repository_root, working_dir, platform, config):
    release_root = config['mounts'].get('release_root')

    branch = config['branch']
    revision = config['revision']

    assert release_root, "config['release_root'] must be set!"

    release_path = f'{release_root}/dummy'

    if not os.path.isdir(release_path): os.makedirs(release_path)

    with open(f'{release_path}/{branch}-{revision}.txt', 'w') as f: f.write('dummy release file\n')

    return {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : f'Config release root set to: {release_root}'}



def run_python_test(repository_root, working_dir, platform, config):

    import zlib, ssl

    python_environment = local_python_install(platform, config)

    if platform['python'][0] == '2': pass
    else:

        if platform['os'] == 'mac' and int( platform['python'].split('.')[1] ) > 6 :
            # SSL certificate test
            import urllib.request; urllib.request.urlopen('https://benchmark.graylab.jhu.edu')

        ves = [
            setup_persistent_python_virtual_environment(python_environment, packages='colr dice xdice pdp11games'),
            setup_python_virtual_environment(working_dir, python_environment, packages='colr dice xdice pdp11games'),
        ]

        for ve in ves:
            commands = [
                'import colr, dice, xdice, pdp11games',
            ]

            if platform['os'] == 'mac' and int( platform['python'].split('.')[1] ) > 6 :
                # SSL certificate test
                commands.append('import urllib.request; urllib.request.urlopen("https://benchmark.graylab.jhu.edu/queue")')

            for command in commands:
                execute('Testing local Python virtual enviroment...', f"{ve.activate} && {ve.python} -c '{command}'")
                execute('Testing local Python virtual enviroment...', f"{ve.activate} && python -c '{command}'")



    return {_StateKey_ : _S_passed_,  _ResultsKey_ : {},  _LogKey_ : f'Done!'}



def compare(test, results, files_path, previous_results, previous_files_path):
    """
    Compare the results of two tests run (new vs. previous) for regression test
    Take two dict and two paths
    Must return standard dict with results

    :param test: str
    :param results: dict
    :param files_path: str
    :param previous_results: dict
    :param previous_files_path: str
    :rtype: dict
    """
    ignore_files = []

    results = dict(tests={}, summary=dict(total=0, failed=0, failed_tests=[]))  # , config={}

    if previous_files_path:
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                exclude = ''.join([' --exclude="{}"'.format(f) for f in ignore_files] ) + ' --exclude="*.ignore"'
                res, brief_diff = execute('Comparing {}...'.format(test), 'diff -rq {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                res, full_diff  = execute('Comparing {}...'.format(test), 'diff -r  {exclude} {0}/{test} {1}/{test}'.format(previous_files_path, files_path, test=test, exclude=exclude), return_='tuple')
                diff = 'Brief Diff:\n' + brief_diff + ( ('\n\nFull Diff:\n' + full_diff[:1024*1024*1]) if full_diff != brief_diff else '' )

                state = _S_failed_ if res else _S_passed_
                results['tests'][test] = {_StateKey_: state, _LogKey_: diff if state != _S_passed_ else ''}

                results['summary']['total'] += 1
                if res: results['summary']['failed'] += 1; results['summary']['failed_tests'].append(test)

    else: # no previous tests case, returning 'passed' for all sub_tests
        for test in os.listdir(files_path):
            if os.path.isdir(files_path + '/' + test):
                results['tests'][test] = {_StateKey_: _S_passed_, _LogKey_: 'First run, no previous results available. Skipping comparison...\n'}
                results['summary']['total'] += 1

    for test in os.listdir(files_path):
        if os.path.isdir(files_path + '/' + test):
            if os.path.isfile(files_path+'/'+test+'/.test_did_not_run.log')  or  os.path.isfile(files_path+'/'+test+'/.test_got_timeout_kill.log'):
                results['tests'][test][_StateKey_] = _S_script_failed_
                results['tests'][test][_LogKey_] += '\nCompare(...): Marking as "Script failed" due to presense of .test_did_not_run.log or .test_got_timeout_kill.log file!\n'
                if test not in results['summary']['failed_tests']:
                    results['summary']['failed'] += 1
                    results['summary']['failed_tests'].append(test)

    state = _S_failed_ if results['summary']['failed'] else _S_passed_

    return {_StateKey_: state, _LogKey_: 'Comparison dummy log...', _ResultsKey_: results}


def run(test, repository_root, working_dir, platform, config, hpc_driver=None, verbose=False, debug=False):
    if   test == 'state':      return run_state_test      (repository_root, working_dir, platform, config)
    elif test == 'regression': return run_regression_test (repository_root, working_dir, platform, config)
    elif test == 'subtests':   return run_subtests_test   (repository_root, working_dir, platform, config)
    elif test == 'release':    return run_release_test    (repository_root, working_dir, platform, config)
    elif test == 'python':     return run_python_test     (repository_root, working_dir, platform, config)

    else: raise BenchmarkError(f'Dummy test script does not support run with test={test!r}!')
