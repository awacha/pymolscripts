#!/usb/bin/env python
import os

from setuptools import setup, find_packages

setup(name='pymolscripts', author='Andras Wacha',
      author_email='awacha@gmail.com', url='http://github.com/awacha/pymolscripts',
      description='Various scripts for pymol',
      packages=find_packages(),
      install_requires=['networkx', 'matplotlib'],
      setup_requires=['setuptools_scm'],
      use_scm_version=True,
      license="BSD 3-clause",
      zip_safe=False,
      )
