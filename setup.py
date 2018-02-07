#!/usr/bin/env python

from distutils.core import setup

setup(
    name='GenomeLocator',
    version='0.1',
    description='Genome Locator',
    author='Jemma Nelson',
    author_email='nelsonjs@altius.org',
    packages=['genome_locator'],
    scripts=['genome_locator/build_genome_index.py'],
    install_requires=[
        'twobitreader',
        'numpy==1.13.3',
        'h5py',
    ],

)
