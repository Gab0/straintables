#!/bin/python

import os


def readGenomeFolder(genomeDirectory="genomes"):
    if os.path.isdir(genomeDirectory):
        genomes = os.listdir(genomeDirectory)
        genomeFilePaths = [os.path.join(genomeDirectory, genomeFile)
                           for genomeFile in genomes
                           if genomeFile.endswith(('.fna', '.fasta'))]

        return genomeFilePaths
    else:
        return []
