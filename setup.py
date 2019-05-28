#!/usr/bin/env python3
#
# Copyright (C) 2019
#

__author__    = 'Mobidic'
__authors__   = ['Henri Pegeot','David Baux','Kevin Yauy','Charles Van Goethem','Thomas Guignard','Olivier Ardouin']
__copyright__ = 'Copyright (C) 2019'
__license__   = 'GNU General Public License v3 (GPLv3)'
__version__   = '2.0a1'
__email__     = 'henri.pegeot@ext.inserm.fr'
__status__    = 'dev'

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="mobidic-MobiCNV",
    version=__version__,
    author=__author__,
    author_email=__email__,
    description="MobiCNV: CNV analysis based on the depth of coverage of Illumina data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mobidic/MobiCNV",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Intended Audience :: Science/Research"
    ],
    install_requires=[
        'pyvcf==0.6.8',
        'numpy==1.16.3',
        'XlsxWriter==1.1.8'
    ],
    entry_points={
        "console_scripts": [
            "mobicnv_main=mobidic_mobicnv:main"
        ],
    },
    scripts = ['scripts/MobiCNV'],
    project_urls={  # Optional
        'Bug Reports': 'https://github.com/mobidic/MobiCNV/issues',
        'Source': 'https://github.com/mobidic/MobiCNV',
},
)
