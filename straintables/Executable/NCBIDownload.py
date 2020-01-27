#!/bin/python

"""

Downloads assemblies found on 'assembly' entrez database,
 they are fetched from 'nucleotide' database.

"""

from Bio import Entrez, SeqIO
import gzip
import json
import os
import re

import random

import shutil

import ftplib
import argparse
import time

from straintables.Database import StrainNames
import socket
# is this legal in the terms of the LAW?
Entrez.email = 'researcher_%i@one-time-use.cn' % random.randrange(0, 1000)

Description = """

straintables' genome and annotation download utility.
Will download genomes and annotations for selected organism
directly from NCBI's FTP servers,
placing them in the correct directories.

It is recommended to select a suitable annotation,
which is defined by the '--strain' argument.

"""


class FTPConnection():
    def __init__(self, ftpServerAddress):
        self.ftpServerAddress = ftpServerAddress

        self.cdir = None

        # ftp connection settings;
        self.timeout = 180
        self.retries = 5
        self.launch()

    def launch(self):
        self.f = ftplib.FTP(self.ftpServerAddress, timeout=self.timeout)
        self.f.login()

        if self.cdir:
            self.cd(self.cdir)

    def cd(self, ftpDirectory):
        self.cdir = ftpDirectory
        self.execute(self.f.cwd, ftpDirectory)

    def listdir(self):
        return self.execute(self.f.nlst)

    def downloadFile(self,
                     remoteFileName,
                     localfilepath):

        with open(localfilepath, 'wb') as localFile:

            downloaded = self.f.retrbinary(
                "RETR %s" % (remoteFileName),
                localFile.write,
                blocksize=32 * 1024)

        return downloaded

    def execute(self, operation, *args):
        for k in range(self.retries):
            if k == self.retries - 1:
                # print("Failure... re-running this script should solve.")
                # exit(1)
                return False
            try:
                result = operation(*args)
                return result

            except (EOFError, TimeoutError,
                    socket.timeout, ftplib.error_temp) as e:
                print("Failure to execute %s on try #%i retries." %
                      (operation.__name__,
                       k + 1))
                print(e)
                try:
                    self.f.quit()
                except Exception:
                    time.sleep(5)
                    self.execute(self.launch)
                    time.sleep(3)
                    print("Retrying....")

            except BrokenPipeError as e:
                print("FTP connection was close due to unknown reasons.")
                print("Operation aborted... please try again.")
                print("\n(%s)\n" % e)
                return False

        return result


class DownloadQuery():
    def __init__(self, IDs, downloadDirectory, FileTypes):
        self.IDs = IDs
        self.downloadDirectory = downloadDirectory
        self.FileTypes = FileTypes

    def execute(self):
        DownloadSuccess = False

        for d, ID in enumerate(self.IDs):
            print("Downloading %i of %i.\n" % (d + 1, len(self.IDs)))
            if not ID:
                DownloadSuccess = False
                continue
            FileSuccess = downloadAssembly(
                ID,
                downloadDirectory=self.downloadDirectory,
                onlyCompleteGenome=True,
                wantedFileTypes=self.FileTypes
            )
            if FileSuccess:
                DownloadSuccess = True

        return DownloadSuccess


def debug(message):
    print('Debug ---> %s' % message)


def findAssemblyList(Organism, Strain=None, retmax=100):
    query = '("%s"[Organism] AND "latest"[Filter]' % Organism
    if Strain:
        query += ' AND "%s"[Strain]' % Strain
    result = Entrez.esearch(db='assembly', term=query, retmax=retmax)

    P = Entrez.read(result)

    # print(json.dumps(P, indent=2))

    return P['IdList']


def downloadAssembly(ID,
                     downloadDirectory='',
                     showSummaries=False,
                     onlyCompleteGenome=True,
                     wantedFileTypes=[],
                     Verbose=False):
    print(ID)
    Summary = Entrez.esummary(db="assembly", id=ID, report='full')
    Summary = Entrez.read(Summary, validate=False)
    relevantSummary = Summary['DocumentSummarySet']['DocumentSummary'][0]

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
        print("Genome Representation is %s" % genomeRepresentation)
        if genomeRepresentation not in ["full", "false"]:
            debug("Skipping partial genome.")
            return 0

    ftpPath = relevantSummary['FtpPath_GenBank']

    debug(ftpPath)

    ftpDirectory = re.findall("genome.*", ftpPath)
    if not ftpDirectory:
        debug("No ftp directory associated. Skipping...")
        return False
    else:
        ftpDirectory = ftpDirectory[0]

    ftpServerAddress = re.findall(r"ftp\..*gov", ftpPath)[0]
    debug(ftpServerAddress)

    # -- instantiate FTP connection;
    Connection = FTPConnection(ftpServerAddress)
    Connection.cd(ftpDirectory)

    # FETCH FTP DIRECTORY FILE LIST;
    remoteFiles = Connection.listdir()
    if not remoteFiles:
        return False

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

    # download only the file with shortest name
    # among matched files at FTP directory;
    remoteFileNames = sorted(remoteFileNames,
                             key=lambda x: len(x))[0: 1]

    showRemoteFiles = True
    if not remoteFileNames or showRemoteFiles:
        if remoteFileNames:
            print("Found %s." % ' '.join(remoteFileNames))
        else:
            print("Remote file key %s not found on FTP server." %
                  ' '.join(wantedFileTypes))
            return 0

        # -- Show list of remote files for accession;
        maxLength = max([len(f) for f in remoteFiles])
        print()
        Header = "Remote files: "
        Header += "=" * (maxLength - len(Header))
        print('\n'.join(remoteFiles))
        print()
        print("=" * maxLength)
        print()

    DownloadSuccess = False
    for remoteFileName in remoteFileNames:
        temporaryFileName = "downloading_" + remoteFileName

        localFinalFilepath = os.path.join(downloadDirectory, remoteFileName)

        if os.path.isfile(localFinalFilepath):
            print("%s already existis! Skipping..." % localFinalFilepath)

        localDownloadingFilepath = os.path.join(downloadDirectory,
                                                temporaryFileName)
        debug("Writing file to %s" % localFinalFilepath)

        debug("Downloading %s" % remoteFileName)

        downloaded = Connection.execute(
            Connection.downloadFile, remoteFileName, localDownloadingFilepath)

        debug(remoteFileName)
        debug("download result: %s" % downloaded)

        if not downloaded:
            return False

        shutil.move(localDownloadingFilepath, localFinalFilepath)
        DownloadSuccess = True

        localFinalFilepath = DecompressFile(localFinalFilepath)

    return DownloadSuccess


def DecompressFile(filepath):
    # GUNZIP FILE - decompress;
    if filepath.endswith(".gz"):
        with gzip.open(filepath, 'rb') as gf:
            file_content = gf.read()
            decompressedPath = os.path.splitext(filepath)[0]
            with open(decompressedPath, 'wb') as f:
                f.write(file_content)
            debug("Decompressed to %s." % decompressedPath)
            os.remove(filepath)
            filepath = decompressedPath
            print('\n')

    return filepath


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
    allGenomes = [
        g for g in allGenomes
        if any([g.endswith(extension) for extension in [".fna", ".fasta"]])
    ]

    def orderByVersion(g):
        for version in range(10, 1, -1):
            if "v%i_" in g.lower():
                return version
        return 1

    allGenomes = sorted(allGenomes, key=orderByVersion)
    for genome in allGenomes:
        localfilepath = os.path.join(dirpath, genome)
        genome = list(SeqIO.parse(localfilepath, format='fasta'))

        if not genome:
            print("Error parsing genome!")
            print("Filename: %s" % localfilepath)
            print(genome)
            continue

        genomeDescription = genome[0].description
        strainName = StrainNames.fetchStrainName(genomeDescription, organism)

        if strainName is not None:
            # A '/' inside the strain name (and output filename) is fatal.
            strainName = strainName.replace("/", "-")
            newFilePath = os.path.join(os.path.dirname(localfilepath),
                                       strainName + '.fna')
            shutil.move(localfilepath, newFilePath)
            print("Moving genome %s file to %s" % (localfilepath, newFilePath))


def parse_arguments():
    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument("--nogenome",
                        dest='downloadGenomes',
                        action='store_false',
                        default=True,
                        help="Skip genome downloads.")

    parser.add_argument("--noannotation",
                        dest='downloadAnnotations',
                        action='store_false',
                        default=True,
                        help="Skip annotation downloads.")

    parser.add_argument("--organism",
                        dest="queryOrganism")

    parser.add_argument("--strain",
                        dest="annotationStrain")

    parser.add_argument("--nbgenomes",
                        dest="genomeSearchMaxResults",
                        type=int,
                        default=20)

    parser.add_argument("--nbannotations",
                        dest="annotationSearchMaxResults",
                        type=int,
                        default=1)

    parser.add_argument("--rename",
                        dest="OnlyRename",
                        help="Rename genome/annotations file and exit.")

    parser.add_argument("--dir",
                        dest="WorkingDirectory",
                        help="Directory to where genome and annotation " +
                        "folders should go",
                        default=".")

    return parser.parse_args()


def Execute(options):

    # -- Define output directories;

    # -- Make sure output directories exist;
    requiredDirectories = ["genomes", "annotations"]
    requiredDirectoriesPath = [
        os.path.join(options.WorkingDirectory, Name)
        for Name in requiredDirectories
        ]

    if options.queryOrganism is None:
        print("Please provide an organism name for the search query.")
        exit(1)

    if not os.path.isdir(options.WorkingDirectory):
        print("Fatal: Selected directory at %s does not exist." %
              options.WorkingDirectory)
        exit(1)

    for Dir in requiredDirectoriesPath:
        if not os.path.isdir(Dir):
            os.mkdir(Dir)

    GenomeDirectory, AnnotationDirectory = requiredDirectoriesPath

    # -- Search Assemblies for Organism;
    AssemblyIDs = findAssemblyList(options.queryOrganism,
                                   retmax=options.genomeSearchMaxResults)

    # -- Abort if no organism is found.
    if not AssemblyIDs:
        print("The search for '%s' returned 0 results... unable to proceed." %
              options.queryOrganism)
        exit(1)

    # -- DOWNLOAD GENOMES;
    dataTypes = []

    # Fetch Genome IDs;
    if options.downloadGenomes:
        dataTypes.append(DownloadQuery(AssemblyIDs,
                                       GenomeDirectory,
                                       ["_genomic.fna"]))

    # Download genomes
    GenomeDownloadSuccess = [query.execute() for query in dataTypes]

    AnnotationDownloadSuccess = []
    # -- DOWNLOAD ANNOTATIONS;
    print("\nFetching Annotations.\n")
    dataTypes = []
    # Fetch Annotation IDs;
    if options.downloadAnnotations:
        if options.annotationStrain:
            # Try to find the user-defined annotation;
            # Fetch user-defined annotation;
            AnnotationIDs = findAssemblyList(
                options.queryOrganism,
                Strain=options.annotationStrain,
                retmax=options.genomeSearchMaxResults
            )

            if AnnotationIDs:
                print("Found Annotation IDs:")
                for i in AnnotationIDs:
                    print(i)
                print()

                dataTypes.append(DownloadQuery(AnnotationIDs,
                                               AnnotationDirectory,
                                               ["_genomic.gbff"]))
                dataTypes.append(DownloadQuery(AnnotationIDs,
                                               GenomeDirectory,
                                               ["_genomic.fna"]))
            else:
                print("User defined annotation strain %s not found!" %
                      options.annotationStrain)
                print("Aborting...")
                exit(1)

        else:
            # Try to download a genome that has a matching annotation;
            print("Annotation not found.")
            dataTypes.append(
                DownloadQuery(
                    [AssemblyIDs[:options.annotationSearchMaxResults]],
                    AnnotationDirectory,
                    ["_genomic.gbff"]
                )
            )

        AnnotationDownloadSuccess = [query.execute() for query in dataTypes]

    # -- REPORT SUCCESS (any?);
    print("\nSummary:")
    if any(GenomeDownloadSuccess):
        print("Genome files downloaded.")

        # properly name genome files;
        renameGenomeFiles(GenomeDirectory, options.queryOrganism)

    if any(AnnotationDownloadSuccess):
        print("Annotation files downloaded.")

    return True


def main():
    options = parse_arguments()

    while True:
        if Execute(options):
            break
        time.sleep(20)


if __name__ == "__main__":
    main()
