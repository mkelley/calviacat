#!/usr/bin/env python3
from glob import glob
from distutils.core import setup

setup(name='calviacat',
      version='1.1.0',
      description='Calibrate star photometry by comparison to a catalog.',
      author='Michael S. P. Kelley',
      author_email='msk@astro.umd.edu',
      url='https://github.com/mkelley/calviacat',
      packages=['calviacat'],
      install_requires=['astropy', 'numpy', 'requests', 'astroquery'],
      license='MIT'
      )
