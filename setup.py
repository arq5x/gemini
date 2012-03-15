import ez_setup
ez_setup.use_setuptools()

import glob
import os
import sys
from setuptools import setup
from distutils.extension import Extension

if 'setuptools.extension' in sys.modules:
    m = sys.modules['setuptools.extension']
    m.Extension.__dict__ = m._Extension.__dict__

version_py = os.path.join(os.path.dirname(__file__), 'pop', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``pop`` is a database framework for exploring genetic variation'
"""

setup(
        name="pop",
        version=version,
        install_requires=['numpy>=1.6.0', 'pyparsing>=1.5.6', 'pysam>=0.6', 'pyvcf>=0.4.2'],
        requires = ['python (>=2.5, <3.0)'],
        packages=['pop',
                  'pop.scripts'],
        author="Aaron Quinlan and Uma Paila",
        description='A database framework for exploring genetic variation',
        long_description=long_description,
        url="none",
        package_dir = {'pop': "pop"},
        scripts = ['pop/scripts/pop'],
        author_email="arq5x@virginia.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
