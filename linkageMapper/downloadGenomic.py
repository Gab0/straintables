#!/bin/python
"""

Downloads assemblies found on 'assembly' entrez database, fetched from 'nucleotide' database.

"""
from Bio import Entrez

import gzip
import json
import os
import re

from ftplib import FTP

Entrez.email = 'probiotics@gut.com'


def debug(message):
    print('Debug ---> %s' % message)


def findAssemblyList(Organism, Strain=None):
    query = '(%s[Organism]' % Organism
    if Strain:
        query += " AND %s[Strain]" % Strain
    result = Entrez.esearch(db='assembly', term=query)

    P = Entrez.read(result)

    print(json.dumps(P, indent=2))
    return P['IdList']


def downloadAssembly(ID, downloadDirectory='',
                     showSummaries=False,
                     onlyCompleteGenome=False,
                     wantedFileTypes=[]):
    print(ID)
    Summary = Entrez.esummary(db="assembly", id=ID, report='full')
    Summary = Entrez.read(Summary, validate=False)
    relevantSummary = Summary['DocumentSummarySet']['DocumentSummary'][0]

    # return Summary
    nucl_id = relevantSummary['AssemblyAccession']

    # print(json.dumps(relevantSummary, indent=2))

    if showSummaries:
        debug("Assembly Summary %s" % json.dumps(Summary, indent=4))
    assemblyStatus = relevantSummary['AssemblyStatus']
    genomeRepresentation = relevantSummary["PartialGenomeRepresentation"]

    debug(assemblyStatus)

    entrySpecies = relevantSummary["Organism"]

    strainData = relevantSummary["Biosource"]["InfraspeciesList"]

    entryStrain = strainData[0]["Sub_value"] if strainData else ''
    print("%s %s" % (entrySpecies, entryStrain))

    print()
    print()

    if onlyCompleteGenome:
        if genomeRepresentation != "false":
            debug("Skipping partial genome.")
            return 0
        """
        if assemblyStatus != 'Complete Genome':
            debug("Skipping scaffold entry.")
            return 0
        """
    ftpPath = relevantSummary['FtpPath_GenBank']
    ftpDirectory = re.findall("genome.*", ftpPath)[0]
    ftpServerAddress = re.findall("ftp\..*gov", ftpPath)[0]
    debug(ftpServerAddress)
    f = FTP(ftpServerAddress)
    f.login()
    f.cwd(ftpDirectory)

    # FETCH FTP DIRECTORY FILE LIST;
    remoteFiles = f.nlst()

    if not wantedFileTypes:
        wantedFileTypes = [
            '_feature_table',
            '_genomic.gff',
            '_genomic.gbff',
            '_genomic.fna'
        ]

    remoteFileNames = []
    for File in remoteFiles:
        for wantedFileType in wantedFileTypes:
            if wantedFileType in File:
                remoteFileNames.append(File)

    remoteFileNames = sorted(remoteFileNames,
                             key=lambda x: len(x))[0:]

    showRemoteFiles = True
    if not remoteFileNames or showRemoteFiles:
        if remoteFileNames:
            print("Found %s." % ' '.join(remoteFileNames))
        else:
            print("Remote file key %s not found on FTP server." % ' '.join(wantedFileTypes))
            return 0
        print('\n'.join(remoteFiles))

    for remoteFileName in remoteFileNames:
        localfilepath = os.path.join(downloadDirectory, remoteFileName)
        debug("Writing assembly to %s" % localfilepath)
        localFile = open(localfilepath, 'wb')

        debug("Downloading %s" % remoteFileName)
        downloaded = f.retrbinary("RETR %s" % (remoteFileName),
                                  localFile.write)

        debug(remoteFileName)
        debug("download result: %s" % downloaded)
        localFile.close()

        # GUNZIP FILE;
        if localfilepath.endswith(".gz"):
            with gzip.open(localfilepath, 'rb') as gf:
                file_content = gf.read()
                decompressedPath = os.path.splitext(localfilepath)[0]
                open(decompressedPath, 'wb').write(file_content)
                debug("Decompressed to %s." % decompressedPath)
                os.remove(localfilepath)
                print('\n')


def strainToDatabase(species):
    P = findAssemblyList(species)
    debug("number of entries: %i" % len(P))
    GenomeCount = 0
    ChromosomeCount = 0

    for Entry in P:
        print("*" * 17)
        n = downloadAssembly(Entry, downloadDirectory='../Analysis/')
        if n:
            GenomeCount += 1
            ChromosomeCount += n

    print("Found %i chromosomes on %i genomes for %s" % (
        ChromosomeCount,
        GenomeCount,
        species))


if __name__ == "__main__":
    D = findAssemblyList("Toxoplasma gondii")
    A = findAssemblyList("Toxoplasma gondii ME49")
    for d in D:
        downloadAssembly(d,
                         downloadDirectory="genomes",
                         onlyCompleteGenome=True,
                         wantedFileTypes=["_genomic.fna"])
    for a in A:
        downloadAssembly(a,
                         downloadDirectory="annotations",
                         onlyCompleteGenome=True,
                         wantedFileTypes=["_genomic.gff"])
