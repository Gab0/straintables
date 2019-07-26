#!/bin/python
import os
from setuptools import setup, find_packages
#from distutils.core import setup


entry_points = {
    'console_scripts': [
        "stpline=straintables.Pipeline:main",
        "stview=straintables.MatrixViewer:main",
        "stdownload=straintables.fetchDataNCBI:main",
        "stprimer=straintables.initializePrimerFile:main"
        ]
}

base_folder = os.path.dirname(os.path.realpath(__file__))
requirements = list(open(os.path.join(base_folder, "requirements.txt")).readlines())
setup(
    name='straintables',
    version='0.93',
    description='Genomic similarities per region',
    author='Gabriel Araujo',
    author_email='gabriel_scf@hotmail.com',
    url='https://www.github.com/Gab0/straintables',
    #packages=find_packages(),
    setup_requires=["numpy"],
    install_requires=requirements,
    packages=[
        'straintables',
        'straintables.Viewer',
        'straintables.PrimerEngine',
        'straintables.DrawGraphics',
        'straintables.Database',
        'straintables.skdistance'
    ],
    platforms='any',
    entry_points=entry_points
)
