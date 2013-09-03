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

version_py = os.path.join(os.path.dirname(__file__), 'gemini', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``gemini`` is a database framework for exploring genetic variation'
"""

setup(
        name="gemini",
        version=version,
        install_requires=['numpy>=1.6.0',
                          'pyparsing>=1.5.6,<=1.5.7',
                          'PyYAML >= 3.10',
                          'cyvcf>=0.1.8',
                          'PyYAML >= 3.10',
                          'python-graph-core >= 1.8.2',
                          'python-graph-dot >= 1.8.2',
                          'bottle >= 0.11.6',
                          'ipython-cluster-helper >= 0.1.9',
                          'bx-python >= 0.7.1',
                          'pandas >= 0.11.0',
                          'scipy >= 0.12.0'],
        dependency_links = ['http://github.com/arq5x/cyvcf/tarball/master#egg=cyvcf-0.1.5'],
        requires = ['python (>=2.5, <3.0)'],
        packages=['gemini',
                  'gemini.scripts',
                  'gemini.data'],
        author="Aaron Quinlan and Uma Paila",
        description='A database framework for exploring genetic variation',
        long_description=long_description,
        url="none",
        package_dir = {'gemini': "gemini"},
        package_data = {'gemini': ['data/gemini.conf']},
        zip_safe = False,
        scripts = ['gemini/scripts/gemini'],
        author_email="arq5x@virginia.edu",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
