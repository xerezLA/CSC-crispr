__author__ = 'Alexendar Perez'

#####################
#                   #
#   Introduction    #
#                   #
#####################

"""CRISPR Specificity Correction Setup Script

National Cancer Institute, National Institutes of Health, United States of America
Developer: Alexendar R. Perez M.D., Ph.D
Primary Investigator: Joana A. Vidigal Ph.D
Laboratory: Vidigal Laboratory, 2019

"""

#################
#               #
#   Libraries   #
#               #
#################

from setuptools import setup, find_packages, Extension
from codecs import open
from os import path

#####################
#                   #
#   Installation    #
#                   #
#####################

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
#'scikit-learn==0.16.1'
setup(
    name = 'CSC_crispr',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['scikit-learn','sklearn-contrib-py-earth==0.1.0','python-dateutil>=2.5.0'],
    dependency_links = ['https://github.com/scikit-learn-contrib/py-earth.git','https://github.com/scikit-learn-contrib/py-earth.git'],
    version = '0.1',
    author= 'Alexendar Perez',
    author_email= 'Alexendar.Perez@ucsf.edu',
    description = 'Computational adjustment for off-targeting effects of NGG Cas9 gRNAs',
    package_data = {'csc_v2' : ['screen_models/Hamming/*.pl','screen_models/examples/*.csv']},
    entry_points={'console_scripts': ['csc_process = csc_v2.csc_lite:main',],})