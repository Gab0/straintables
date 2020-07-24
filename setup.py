#!/bin/python
import os
from setuptools import setup
#from distutils.core import setup


Executables = [
    "stgenomepline=straintables.Executable.GenomePipeline:main",
    "stview=straintables.Executable.WebViewer:main",
    "stdownload=straintables.Executable.NCBIDownload:main",
    "stgenregion=straintables.Executable.initializePrimerFile:main",
    "stprotein=straintables.Executable.Protein:main",
    "stfastapline=straintables.Executable.FastaPipeline:main",
    "stblast=straintables.Executable.BlastAmplicons:main"
]

Tools = []

entry_points = {
    'console_scripts':  Executables + Tools
}

base_folder = os.path.dirname(os.path.realpath(__file__))
requirements = list(
    open(os.path.join(base_folder, "requirements.txt")).readlines())

setup(
    name='straintables',
    version='1.53',
    description='Build & Compare dissimilarity matrices for genomic regions',
    author='Gabriel Araujo',
    author_email='gabriel_scf@hotmail.com',
    url='https://www.github.com/Gab0/straintables',
    # packages=find_packages(),
    setup_requires=["numpy"],
    install_requires=requirements,
    packages=[
        'straintables',
        'straintables.Executable',
        'straintables.Viewer',
        'straintables.PrimerEngine',
        'straintables.DrawGraphics',
        'straintables.Database',
        'straintables.skdistance',
        'straintables.DistanceEngine'
    ],
    package_data={'': [
        'Viewer/WebComponents/' + f
        for f in [
                'MainView.html',
                'style.css',
                'CustomPlotView.html',
                'CustomPlotBuild.html',
                'Menu.html',
                'Footer.html',
                'PlotOptions.html',
                'dropdown_logic.js'
        ]
    ]},
    include_package_data=True,
    platforms='any',
    entry_points=entry_points
)
