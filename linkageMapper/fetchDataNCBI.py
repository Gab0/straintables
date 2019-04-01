#!/bin/python
"""

Downloads assemblies found on 'assembly' entrez database,
 they are fetched from 'nucleotide' database.

"""
from Bio import Entrez, SeqIO, Seq, Alphabet
import gzip
import json
import os
import re

import random

import shutil

from ftplib import FTP
import optparse

# is this legal in the terms of the LAW?
Entrez.email = 'researcher_%i@one-time-use.cn' % random.randrange(0, 1000)


def debug(message):
    print('Debug ---> %s' % message)


def findAssemblyList(Organism, Strain=None, retmax=100):
    query = '(%s[Organism]' % Organism
    if Strain:
        query += " AND %s[Strain]" % Strain
    result = Entrez.esearch(db='assembly', term=query, retmax=100)

    P = Entrez.read(result)

    print(json.dumps(P, indent=2))
    return P['IdList']


def downloadAssembly(ID,
                     downloadDirectory='',
                     showSummaries=False,
                     onlyCompleteGenome=False,
                     wantedFileTypes=[],
                     Verbose=False):
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

    debug(ftpPath)

    ftpDirectory = re.findall("genome.*", ftpPath)
    if not ftpDirectory:
        debug("No ftp directory associated. Skipping...")
        return False
    else:
        ftpDirectory = ftpDirectory[0]

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

    remoteFileExtension = None
    remoteFileNames = []
    for File in remoteFiles:
        for wantedFileType in wantedFileTypes:
            if wantedFileType in File:
                remoteFileNames.append(File)
                remoteFileExtension = "." + wantedFileType.split(".")[-1]

    # download only the file with shortest name among matched files at FTP directory;
    remoteFileNames = sorted(remoteFileNames,
                             key=lambda x: len(x))[0: 1]

    showRemoteFiles = True
    if not remoteFileNames or showRemoteFiles:
        if remoteFileNames:
            print("Found %s." % ' '.join(remoteFileNames))
        else:
            print("Remote file key %s not found on FTP server." % ' '.join(wantedFileTypes))
            return 0

        print("\nRemote files:")
        print('\n'.join(remoteFiles))
        print()

    for remoteFileName in remoteFileNames:
        localFileName = remoteFileName

        localfilepath = os.path.join(downloadDirectory, localFileName)
        debug("Writing assembly to %s" % localfilepath)
        localFile = open(localfilepath, 'wb')

        debug("Downloading %s" % remoteFileName)
        downloaded = f.retrbinary("RETR %s" % (remoteFileName),
                                  localFile.write)

        debug(remoteFileName)
        debug("download result: %s" % downloaded)
        localFile.close()

        # GUNZIP FILE - decompress;
        if localfilepath.endswith(".gz"):
            with gzip.open(localfilepath, 'rb') as gf:
                file_content = gf.read()
                decompressedPath = os.path.splitext(localfilepath)[0]
                open(decompressedPath, 'wb').write(file_content)
                debug("Decompressed to %s." % decompressedPath)
                os.remove(localfilepath)
                localfilepath = decompressedPath
                print('\n')

        # ANNOTATION FILE - split by scaffold (unused);
        if False and localfilepath.endswith(".gbff"):
            with open(localfilepath) as af:
                data = list(SeqIO.parse(af, 'gb'))
                for chromosome in data:
                    seq_len = len(str(chromosome.seq))
                    print(chromosome.description)

                    label = re.findall("chromosome (\w+),", chromosome.description)
                    if label:
                        label = label[0]
                        outputFile = "chromosome_%s" % label

                        outputFilePath = os.path.join(downloadDirectory, outputFile)
                        print("Creating %s" % outputFilePath)

                        # erase chromosome sequence;
                        chromosome.seq = Seq.Seq("", Alphabet.generic_dna)
                        SeqIO.write(chromosome, outputFilePath, 'gb')

            os.remove(localfilepath)


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


def renameGenomeFiles(dirpath, organism):
    allGenomes = os.listdir(dirpath)

    def orderByVersion(g):
        for version in range(10, 1, -1):
            if "v%i_" in g.lower():
                return version
        return 1

    allGenomes = sorted(allGenomes, key=orderByVersion)
    for genome in allGenomes:
        localfilepath = os.path.join(dirpath, genome)
        genome = list(SeqIO.parse(localfilepath, format='fasta'))

        genomeDescription = genome[0].description

        strainName = None

        # RENAMING CRITERIA 1:
        def renamingCriteria1():
            d = re.findall("strain (\w+)", genomeDescription)
            if d:
                print("C1")
                return d[0]

        # RENAMING CRITERIA 2:
        def renamingCriteria2():
            e = re.findall("%s ([\w-\d]+)" % (organism), genomeDescription, flags=re.IGNORECASE)
            if e:
                return e[0]

        # RENAMING CRITERIA 3:
        def renamingCriteria3():
            words = genomeDescription.split(" ")[1:]
            uppercaseWords = [w.isupper() for w in words]
            if sum(uppercaseWords) == 1:
                return words[uppercaseWords.index(True)]

        allCriteria = [
            renamingCriteria2,
            renamingCriteria1,
            renamingCriteria3
        ]

        for renamingCriteria in allCriteria:
            strainName = renamingCriteria()
            if strainName is not None:
                newFilePath = os.path.join(os.path.dirname(localfilepath), strainName + '.fna')
                shutil.move(localfilepath, newFilePath)
                print("Moving genome %s file to %s" % (localfilepath, newFilePath))
                break


def parse_options():
    parser = optparse.OptionParser()
    parser.add_option("--nogenome", dest='downloadGenomes',
                      action='store_false',
                      default=True,
                      help="Skip genome downloads.")

    parser.add_option("--noannotation",
                      dest='downloadAnnotations',
                      action='store_false',
                      default=True,
                      help="Skip annotation downloads.")

    parser.add_option("--organism",
                      dest="queryOrganism",
                      default="Toxoplasma gondii")

    parser.add_option("--strain",
                      dest="annotationStrain",
                      default="ME49")

    parser.add_option("--maxgenomecount",
                      dest="genomeSearchMaxResults",
                      type=int,
                      default=20)

    options, args = parser.parse_args()

    return options


def main():
    options = parse_options()
    dataTypes = []

    # Fetch Genome IDs;
    if options.downloadGenomes:
        D = findAssemblyList(options.queryOrganism, retmax=options.genomeSearchMaxResults)
        dataTypes.append((D, "genomes", ["_genomic.fna"]))

    # -- Make sure output directories exist;
    requiredDirectories = ["genomes", "annotations"]
    for Dir in requiredDirectories:
        if not os.path.isdir(Dir):
            os.mkdir(Dir)

    # Fetch Annotation IDs;
    if options.downloadAnnotations:
        A = findAssemblyList("%s %s" % (options.queryOrganism, options.annotationStrain))
        dataTypes.append((A, "annotations", ["_genomic.gbff"]))

    for (IDS, typeName, fileExtensions) in dataTypes:
        for d, ID in enumerate(IDS):
            print("Downloading %i of %i.\n" % (d + 1, len(IDS)))

            downloadAssembly(ID,
                             downloadDirectory=typeName,
                             onlyCompleteGenome=True,
                             wantedFileTypes=fileExtensions)

    # properly name genome files;
    renameGenomeFiles("genomes", options.queryOrganism)

    print()
    print("Sucess:")
    print("Annotation and genome files downloaded.")


if __name__ == "__main__":
    main()
