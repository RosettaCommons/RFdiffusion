from setuptools import setup, find_packages

setup(name='rfdiffusion',
      version='1.1.0',
      description='RFdiffusion is an open source method for protein structure generation.',
      author='Rosetta Commons',
      url='https://github.com/RosettaCommons/RFdiffusion',
      scripts=["scripts/run_inference.py"],
      packages=find_packages(),
      install_requires=['torch', 'se3-transformer'])
