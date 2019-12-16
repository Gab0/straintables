#!/bin/python
import os
from setuptools import setup
#from distutils.core import setup


entry_points = {
    'console_scripts': [
        "stpline=straintables.Executable.Pipeline:main",
        "stview=straintables.Executable.WebViewer:main",
        "stdownload=straintables.Executable.fetchDataNCBI:main",
        "stgenregion=straintables.Executable.initializePrimerFile:main",
        "stprotein=straintables.Executable.Protein:main",
        "stfromfasta=straintables.Executable.fromMultifasta:main"
        ]
}

base_folder = os.path.dirname(os.path.realpath(__file__))
requirements = list(
    open(os.path.join(base_folder, "requirements.txt")).readlines())

setup(
    name='straintables',
    version='1.3',
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
    ],
    package_data={'': [
        'Viewer/WebComponents/' + f
        for f in [
                'MainView.html',
                'style.css',
                'CustomPlotView.html',
                'CustomPlotBuild.html',
                'Menu.html',
                'PlotOptions.html',
                'dropdown_logic.js'
        ]
    ]},
    include_package_data=True,
    platforms='any',
    entry_points=entry_points
)
