#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   benchmark.py
## @brief  Run arbitrary Rosetta testing script
## @author Sergey Lyskov

from __future__ import print_function

import os, os.path, sys, shutil, json, platform, re
import codecs

from importlib.machinery import SourceFileLoader

from configparser import ConfigParser, ExtendedInterpolation
import argparse

from tests import *  # execute, Tests states and key names
from hpc_drivers import *


# Calculating value of Platform dict
Platform = {}
if sys.platform.startswith("linux"):
    Platform['os'] = 'ubuntu' if os.path.isfile('/etc/lsb-release') and 'Ubuntu' in open('/etc/lsb-release').read() else 'linux'  # can be linux1, linux2, etc
elif sys.platform == "darwin" :      Platform['os'] = 'mac'
elif sys.platform == "cygwin" :      Platform['os'] = 'cygwin'
elif sys.platform == "win32" :       Platform['os'] = 'windows'
else:                                Platform['os'] = 'unknown'

#Platform['arch'] = platform.architecture()[0][:2]  # PlatformBits
Platform['compiler'] = 'gcc' if Platform['os'] == 'linux' else 'clang'

Platform['python'] = sys.executable


def load_python_source_from_file(module_name, module_path):
    ''' replacment for deprecated imp.load_source
    '''
    return SourceFileLoader(module_name, module_path).load_module()


class Setup(object):
    __slots__ = 'test working_dir platform config compare debug'.split()  # version daemon path_to_previous_test
    def __init__(self, **attrs):
        #self.daemon = True
        for k, v in attrs.items():
            if k in self.__slots__: setattr(self, k, v)


def setup_from_options(options):
    ''' Create Setup object based on user supplied options, config files and auto-detection
    '''
    platform = dict(Platform)

    if options.suffix: options.suffix = '.' + options.suffix

    platform['extras'] = options.extras.split(',') if options.extras else []
    platform['python'] = options.python
    #platform['options'] = json.loads( options.options ) if options.options else {}

    if options.memory: memory = options.memory
    elif platform['os'] in ['linux', 'ubuntu']: memory = int( execute('Getting memory info...', 'free -m', terminate_on_failure=False, silent=True, silence_output_on_errors=True, return_='output').split('\n')[1].split()[1]) // 1024
    elif platform['os'] == 'mac':   memory = int( execute('Getting memory info...', 'sysctl -a | grep hw.memsize', terminate_on_failure=False, silent=True, silence_output_on_errors=True, return_='output').split()[1]) // 1024 // 1024 // 1024

    platform['compiler'] = options.compiler

    if os.path.isfile(options.config):
        with open(options.config) as f:
            if '%(here)s' in f.read():
                print(f"\n\n>>> ERROR file `{options.config}` seems to be in outdated format! Please use benchmark.template.ini to update it.")
                sys.exit(1)

        user_config = ConfigParser(
            dict(
                _here_ = os.path.abspath('./'),
                _user_home_ = os.environ['HOME']
            ),
            interpolation = ExtendedInterpolation()
        )

        with open(options.config) as f: user_config.readfp(f)

    else:
        print(f"\n\n>>> Config file `{options.config}` not found. You may want to manually copy `benchmark.ini.template` to `{options.config}` and edit the settings\n\n")
        user_config = ConfigParser()
        user_config.set('main', 'cpu_count',  '1')
        user_config.set('main', 'hpc_driver', 'MultiCore')
        user_config.set('main', 'branch',     'unknown')
        user_config.set('main', 'revision',   '42')
        user_config.set('main', 'user_name',  'Jane Roe')
        user_config.set('main', 'user_email', 'jane.roe@university.edu')
        user_config.add_section('main')

    if options.jobs: user_config.set('main', 'cpu_count', str(options.jobs) )
    user_config.set('main', 'memory',    str(memory) )

    if options.mount:
        for m in options.mount:
            key, _, path = m.partition(':')
            user_config.set('mount', key, path)

    #config = Config.items('config')
    #for section in config.sections(): print('Config section: ', section, dict(config.items(section)))
    #config = { section: dict(Config.items(section)) for section in Config.sections() }

    config = { k : d for k, d in user_config['main'].items() if k not in user_config[user_config.default_section] }
    config['mounts'] = { k : d for k, d in user_config['mount'].items() if k not in user_config[user_config.default_section] }

    #print(json.dumps(config, sort_keys=True, indent=2)); sys.exit(1)

    #config.update( config.pop('config').items() )

    config = dict(config,
                  cpu_count = user_config.getint('main', 'cpu_count'),
                  memory = memory,
                  revision = user_config.getint('main', 'revision'),
                  emulation=True,
    )  # debug=options.debug,

    if 'results_root' not in config: config['results_root'] = os.path.abspath('./results/')

    if 'prefix' in config:
        assert os.path.isabs( config['prefix'] ), f'ERROR: `prefix` path must be absolute! Got: {config["prefix"]}'

    else: config['prefix'] = os.path.abspath( config['results_root'] + '/prefix')

    config['merge_head'] = options.merge_head
    config['merge_base'] = options.merge_base

    if options.skip_compile is not None: config['skip_compile'] = options.skip_compile

    #print(f'Results path: {config["results_root"]}')
    #print('Config:{}, Platform:{}'.format(json.dumps(config, sort_keys=True, indent=2), Platform))

    if options.compare: print('Comparing tests {} with suffixes: {}'.format(options.args, options.compare) )
    else: print('Running tests: {}'.format(options.args) )

    if len(options.args) != 1: print('Error: Single test-name-to-run should be supplied!');  sys.exit(1)
    else:
        test = options.args[0]
        if test.startswith('tests/'): test = test.partition('tests/')[2][:-3]  # removing dir prefix and .py suffix

    if options.compare:
        compare = options.compare[0], options.compare[1]  # (this test suffix, previous test suffix)
        working_dir = os.path.abspath( config['results_root'] + f'/{platform["os"]}.{test}' )  # will be a root dir with sub-dirs (options.compare[0], options.compare[1])
    else:
        compare = None
        working_dir = os.path.abspath( config['results_root'] + f'/{platform["os"]}.{test}{options.suffix}' )


    if os.path.isdir(working_dir): shutil.rmtree(working_dir);  #print('Removing old job dir %s...' % working_dir)  # remove old dir if any
    os.makedirs(working_dir)

    setup = Setup(
        test        = test,
        working_dir = working_dir,
        platform    = platform,
        config      = config,
        compare     = compare,
        debug       = options.debug,
        #daemon      = False,
    )

    setup_as_json = json.dumps( { k : getattr(setup, k) for k in setup.__slots__}, sort_keys=True, indent=2)
    with open(working_dir + '/.setup.json', 'w') as f: f.write(setup_as_json)

    #print(f'Detected hardware platform: {Platform}')
    print(f'Setup: {setup_as_json}')
    return setup


def truncate_log(log):
    _max_log_size_  = 1024*1024*1
    _max_line_size_ = _max_log_size_ // 2

    if len(log) > _max_log_size_:
        new = log
        lines = log.split('\n')

        if len(lines) > 256:
            new_lines = lines[:32] + ['...truncated...'] + lines[-128:]
            new = '\n'.join(new_lines)

        if len(new) > _max_log_size_: # special case for Ninja logs that does not use \n
            lines = re.split(r'[\r\n]*', log)  #t.log.split('\r')
            if len(lines) > 256: new = '\n'.join( lines[:32] + ['...truncated...'] + lines[-128:] )

        if len(new) > _max_log_size_: # going to try to truncate each individual line...
            print(f'Trying to truncate log line-by-line...')
            new = '\n'.join( (
                ( line[:_max_line_size_//3] + '...truncated...' + line[-_max_line_size_//3:] ) if line > _max_line_size_ else line
                for line in new_lines ) )

        if len(new) > _max_log_size_: # fall-back strategy in case all of the above failed...
            print(f'WARNING: could not truncate log line-by-line, falling back to raw truncate...')
            new = 'WARNING: could not truncate test log line-by-line, falling back to raw truncate!\n...truncated...\n' + ( '\n'.join(lines) )[-_max_log_size_+256:]

        print( 'Trunacting test output log: {0}MiB --> {1}MiB'.format(len(log)/1024/1024, len(new)/1024/1024) )

        log = new

    return log

def truncate_results_logs(results):
    results[_LogKey_] = truncate_log( results[_LogKey_] )
    if _ResultsKey_ in results  and  _TestsKey_ in results[_ResultsKey_]:
        tests = results[_ResultsKey_][_TestsKey_]
        for test in tests:
            tests[test][_LogKey_] = truncate_log( tests[test][_LogKey_] )


def find_test_description(test_name, test_script_file_name):
    ''' return content of test-description file if any or None if no description was found
    '''

    def find_description_file(prefix, test_name):
        fname = prefix + test_name + '.md'
        if os.path.isfile(fname): return fname
        return prefix + 'md'

    description_file_name =  find_description_file( test_script_file_name[:-len('command.py')] + 'description.', test_name) if test_script_file_name.endswith('/command.py') else find_description_file(test_script_file_name[:-len('py')], test_name)

    if description_file_name  and  os.path.isfile(description_file_name):
        print(f'Found test suite description in file: {description_file_name!r}')
        with open(description_file_name, encoding='utf-8', errors='backslashreplace') as f: description = f.read()
        return description

    else: return None



def run_test(setup):
    #print(f'{setup!r}')
    suite, rest = setup.test.split('.'), []
    while suite:
        #print( f'suite: {suite}, test: {rest}' )

        file_name = '/'.join( ['tests'] + suite ) + '.py'
        if os.path.isfile(file_name): break

        file_name = '/'.join( ['tests'] + suite ) + '/command.py'
        if os.path.isfile(file_name): break

        rest.insert(0, suite.pop())


    test = '.'.join( suite + rest )
    test_name = '.'.join(rest)

    print( f'Loading test from: {file_name}, suite+test: {test!r}, test: {test_name!r}' )
    #test_suite = imp.load_source('test_suite', file_name)
    test_suite = load_python_source_from_file('test_suite', file_name)

    test_description = find_test_description(test_name, file_name)

    if setup.compare:
        #working_dir_1 = os.path.abspath( config['results_root'] + f'/{Platform["os"]}.{test}.{Options.compare[0]}' )
        working_dir_1 = setup.working_dir + f'/{setup.compare[0]}'

        working_dir_2        = setup.compare[1]  and  ( setup.working_dir + f'/{setup.compare[1]}' )
        res_2_json_file_path = setup.compare[1]  and  f'{working_dir_2}/.execution.results.json'

        with open(working_dir_1 + '/.execution.results.json') as f: res_1 = json.load(f).get(_ResultsKey_)

        if setup.compare[1] and ( not os.path.isfile(res_2_json_file_path) ):
            setup.compare[1] = None
            state_override = _S_failed_
        else:
            state_override = None

        if setup.compare[1] == None: res_2, working_dir_2 = None, None
        else:
            with open(res_2_json_file_path) as f: res_2 = json.load(f).get(_ResultsKey_)

        res = test_suite.compare(test, res_1, working_dir_1, res_2, working_dir_2)

        if state_override:
            log_prefix = \
                f'WARNING: Previous test results does not have `.execution.results.json` file, so comparision with None was performed instead!\n' \
                f'WARNING: Overriding calcualted test state `{res[_StateKey_]}` → `{_S_failed_}`...\n\n'

            res[_LogKey_] = log_prefix + res[_LogKey_]
            res[_StateKey_] = _S_failed_


        # # Caution! Some of the strings in the result object may be unicode. Be robust to unicode in the log messages.
        # with codecs.open(setup.working_dir+'/.comparison.log.txt', 'w', encoding='utf-8', errors='replace') as f: f.write( truncate_log( res[_LogKey_] ) )
        # res[_LogKey_] = truncate_log( res[_LogKey_] )

        # # Caution! Some of the strings in the result object may be unicode. Be robust to unicode in the log messages.
        with codecs.open(setup.working_dir+'/.comparison.log.txt', 'w', encoding='utf-8', errors='replace') as f: f.write(res[_LogKey_])
        truncate_results_logs(res)

        print( 'Comparison finished with output:\n{}'.format( res[_LogKey_] ) )

        with open(setup.working_dir+'/.comparison.results.json', 'w') as f: json.dump(res, f, sort_keys=True, indent=2)

        #print( 'Comparison finished with results:\n{}'.format( json.dumps(res, sort_keys=True, indent=2) ) )
        if 'summary' in res: print('Summary section:\n{}'.format( json.dumps(res['summary'], sort_keys=True, indent=2) ) )

        print( f'Output results of this comparison saved to {working_dir_1}/.comparison.results.json\nComparison log saved into {working_dir_1}/.comparison.log.txt' )


    else:
        working_dir = setup.working_dir  #os.path.abspath( setup.config['results_root'] + f'/{platform["os"]}.{test}{options.suffix}' )

        hpc_driver_name = setup.config['hpc_driver']
        hpc_driver = None if hpc_driver_name in ['', 'none'] else eval(hpc_driver_name + '_HPC_Driver')(working_dir, setup.config, tracer=print, set_daemon_message=lambda x:None)

        api_version = test_suite._api_version_ if hasattr(test_suite, '_api_version_') else ''

        # if api_version < '1.0':
        #     res = test_suite.run(test=test_name, rosetta_dir=os.path.abspath('../..'), working_dir=working_dir, platform=dict(Platform), jobs=Config.cpu_count, verbose=True, debug=Options.debug)
        # else:

        if api_version == '1.0': res = test_suite.run(test=test_name, repository_root=os.path.abspath('./..'), working_dir=working_dir, platform=dict(setup.platform), config=setup.config, hpc_driver=hpc_driver, verbose=True, debug=setup.debug)
        else:
            print(f'Test benchmark api_version={api_version} is not supported!'); sys.exit(1)

        if not isinstance(res, dict): print(f'Test returned result of type {type(res)} while dict-like object was expected, please check that test-script have correct `return` statment! Terminating...'); sys.exit(1)

        # Caution! Some of the strings in the result object may be unicode. Be robust to unicode in the log messages
        with codecs.open(working_dir+'/.execution.log.txt', 'w', encoding='utf-8', errors='replace') as f: f.write( res[_LogKey_] )

        # res[_LogKey_] = truncate_log( res[_LogKey_] )
        truncate_results_logs(res)

        if _DescriptionKey_ not in res: res[_DescriptionKey_] = test_description

        if res[_StateKey_] not in _S_Values_: print( 'Warning!!! Test {} failed with unknow result code: {}'.format(test_name, res[_StateKey_]) )
        else: print( f'Test {test} finished with output:\n{res[_LogKey_]}\n----------------------------------------------------------------\nState: {res[_StateKey_]!r} | ', end='')

        # JSON by default serializes to an ascii-encoded format
        with open(working_dir+'/.execution.results.json', 'w') as f: json.dump(res, f, sort_keys=True, indent=2)

        print( f'Output and full log of this test saved to:\n{working_dir}/.execution.results.json\n{working_dir}/.execution.log.txt' )






def main(args):
    ''' Script to Run arbitrary Rosetta test
    '''
    parser = argparse.ArgumentParser(usage="Main testing script to run tests in the tests directory. "
                                           "Use the --skip-compile to skip the build phase when testing locally. "
                                           "Example Command: /benchmark.py -j2 integration.valgrind")

    parser.add_argument('-j', '--jobs', default=0, type=int, help="Number of processors to use on when building. (default: use value from config file or 1)")

    parser.add_argument('-m', '--memory', default=0, type=int, help="Amount of memory to use (default: use 2Gb per job")

    parser.add_argument('--compiler', default=Platform['compiler'], help="Compiler to use")

    #parser.add_argument('--python', default=('3.9' if Platform['os'] == 'mac' else '3.6'), help="Python interpreter to use")
    parser.add_argument('--python', default=f'{sys.version_info.major}.{sys.version_info.minor}.s', help="Specify version of Python interpreter to use, for example '3.9'. If '.s' added to end of version string then use the same interpreter that was used to start this script. Default: '?.?.s'")

    parser.add_argument("--extras", default='', help="Specify scons extras separated by ',': like --extras=mpi,static" )

    parser.add_argument("--debug", action="store_true", dest="debug", default=False, help="Run specified test in debug mode (not with debug build!) this mean different things and depend on the test. Could be: skip the build phase, skip some of the test phases and so on. [off by default]" )

    parser.add_argument("--suffix", default='', help="Specify ending suffix for test output dir. This is useful when you want to save test results in different dir for later comparison." )

    parser.add_argument("--compare", nargs=2, help="Do not run the tests but instead compare previous results. Use --compare suffix1 suffix2" )

    parser.add_argument("--config", default='benchmark.{os}.ini'.format(os=Platform['os']), action="store", help="Location of .ini file with additional options configuration. Optional.")

    parser.add_argument("--skip-compile", dest='skip_compile', default=None, action="store_true", help="Skip the compilation phase. Assumes the binaries are already compiled locally.")

    #parser.add_argument("--results-root", default=None, action="store", help="Location of `results` dir default is to use `./results`")

    parser.add_argument("--setup", default=None, help="Specify JSON file with setup information. When this option supplied all other config and commandline options is ignored and auto-detection disable. Test, platform info will be gathered from provided JSON file. This option is designed to be used in daemon mode." )

    parser.add_argument("--merge-head", default='HEAD', help="Specify SHA1/branch-name that will be used for `merge-head` value when simulating PR testing" )

    parser.add_argument("--merge-base", default='origin/master', help="Specify SHA1/branch-name that will be used for `merge-base` value when simulating PR testing" )

    parser.add_argument("--mount", action="append", help="Specify one of the mount points, like: --mount release_root:/some/path. This option could be used multiple times if needed" )


    parser.add_argument('args', nargs=argparse.REMAINDER)

    options = parser.parse_args(args=args[1:])

    if any( [a.startswith('-') for a in options.args] ) :
        print( '\nWARNING WARNING WARNING WARNING\n' )
        print( '\tInterpreting', ' '.join(["'"+a+"'" for a in options.args if a.startswith('-')]), 'as test name(s), rather than as option(s).'  )
        print( "\tTry moving it before any test name, if that's not what you want."  )
        print( '\nWARNING WARNING WARNING WARNING\n'  )


    if options.setup:
        with open(options.setup) as f: setup = Setup( **json.load(f) )

    else:
        setup = setup_from_options(options)

    run_test(setup)


if __name__ == "__main__": main(sys.argv)
