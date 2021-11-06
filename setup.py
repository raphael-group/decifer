"""
setup.py

author: Simone Zaccaria
date: 2021-01-15
"""


import setuptools
from setuptools import setup


setuptools.setup(
    name='decifer',
    version='v2.0.3',
    python_requires='>=3.7.*',
    packages=['decifer'],
    package_dir={'': 'src'},
    author='Simone Zaccaria and Gryte Satas and Mohammed El-Kebir',
    author_email='zaccaria@ucl.ac.uk and gsatas7@gmail.com and melkebir@illinois.edu',
    description='DeCiFer is an algorithm that simultaneously selects mutation multiplicities and clusters SNVs by their corresponding descendant cell fractions (DCF).',
    url='https://github.com/raphael-group/decifer',
    install_requires=[
        'numpy>=1.16.1',
        'scipy>=1.2.1',
        'pandas',
        'lemon',
        'seaborn>=0.7.1'
    ],
    include_package_data = True,
    package_data = {
        '' : ['*.txt']
        },
    license='BSD',
    platforms=["Linux", "MacOs", "Windows"],
    classifiers=[
        'Programming Language :: Python :: 3.7',
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX :: Linux",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords=[
        'scientific',
        'sequence analysis',
        'cancer',
        'somatic mutations'
        'DNA',
        'SNVs'],
    entry_points={'console_scripts': ['decifer=decifer.__main__:main',
                                      'fitbeta=decifer.generator:main']}
)
