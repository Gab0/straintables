#!/bin/python
import os
import re
import itertools
import json
import random
import numpy as np
import pandas as pd

import sys
from collections import OrderedDict

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Seq, SeqIO

from optparse import OptionParser

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))


from . import PrimerEngine


class clonalTypeReference():
    def __init__(self):
        self.genotypeData = pd.read_csv("genomes_haplogroups.csv")
        self.rflpGenotypes = pd.read_csv("genotypes.csv")

    def getGenotypeNumber(self, name):
        found = self.genotypeData[self.genotypeData.Genome == name]
        GenotypeNumber = found.iloc[0].ToxoDB

        return GenotypeNumber

    def getRFLPLocus(self, genotypeNumber, referenceLocus):

        found = self.rflpGenotypes[
            self.rflpGenotypes.Genotype == genotypeNumber]
        RFLPLocus = found.iloc[0][referenceLocus]
        return RFLPLocus


def writeFastaFile(outputPath,
                   locusName,
                   locusSequences,
                   clonalReference=None):
    fastaSequences = []

    for genome in locusSequences.keys():
        Name = genome + ".fasta"
        try:
            if clonalReference:
                REF = options.LocusReference
                referenceLocus = REF if REF else options.LocusName

                GN = clonalReference.getGenotypeNumber(Name)
                RFLPLocus = clonalReference.getRFLPLocus(GN, REF)
                Name += "___%s" % RFLPLocus
                print("Loci data found for %s" % Name)
        except Exception as e:
            pass

        sequence = SeqIO.SeqRecord(Seq.Seq(locusSequences[genome]),
                                   id=genome, name=genome, description="")
        fastaSequences.append(sequence)

    with open(outputPath, "w") as output_handle:
        SeqIO.write(fastaSequences, output_handle, "fasta")


# LOAD USER DEFINED PRIMER DATA, with or without header;
def loadPrimerList(filePath):
    lociPrimerList = pd.read_csv(filePath)
    expectedColumns = ["LocusName", "ForwardPrimer", "ReversePrimer"]

    fileColumns = list(lociPrimerList.columns)
    for i in range(len(fileColumns)):
        if "Unnamed" in fileColumns[i]:
            fileColumns[i] = np.nan
    if fileColumns != expectedColumns:

        newFirstRowData = dict([(expected, fileColumns[e])
                                for e, expected in enumerate(expectedColumns)])

        newFirstRow = pd.DataFrame([newFirstRowData], columns=expectedColumns)
        lociPrimerList.columns = expectedColumns

        lociPrimerList = pd.concat([newFirstRow, lociPrimerList],
                                   axis=0, ignore_index=True).reset_index(drop=True)

    return lociPrimerList


def Execute(options):
    # LOAD CLONAL TYPE LOCUS INFORMATION (Su et al.);
    #clonalReference = clonalTypeReference()
    clonalReference = None

    # CHECK DECLARATION OF PRIMER FILE;
    if not options.primerFile:
        print("FATAL: No primer file specified.")
        exit(1)

    PrimerTypes = ["ForwardPrimer", "ReversePrimer"]

    # -- LOAD GENOME FEATURES;
    featureFolderPath = "annotations"
    if os.path.isdir(featureFolderPath):
        genomeFeatureFiles = [
            os.path.join(featureFolderPath, File)
            for File in os.listdir(featureFolderPath)
            if not File.startswith(".")
        ]

        genomeFeatures = [
            list(SeqIO.parse(File, "genbank"))
            for File in genomeFeatureFiles if File.endswith(".gbff")
        ]

    else:
        genomeFeatures = []

    # CHECK GENOME FEATURES FILE EXISTENCE;
    if genomeFeatures:
        genomeFeatures = genomeFeatures[0]
    else:
        print("Fatal: No features found.")

    """
    # FETCH GENOME NAMES;
    genomeFeaturesChromosomes =\
        [
            chromosome.features[0].qualifiers['chromosome'][0].upper()
            for chromosome in genomeFeatures
            if 'chromosome' in chromosome.features[0].qualifiers.keys()
        ]

    genomeFeaturesChromosomes = [
        chrName
        for chrName in genomeFeaturesChromosomes
        if chrName is not "UNKNOWN"
    ]

    # print(genomeFeaturesChromosomes)
    """

    # -- LOAD USER DEFINED PRIMERS;
    lociPrimerList = loadPrimerList(options.primerFile)

    # LOAD GENOMES;
    genomeDirectory = "genomes"
    if os.path.isdir(genomeDirectory):
        genomes = os.listdir(genomeDirectory)
        genomeFilePaths = [os.path.join(genomeDirectory, genomeFile)
                           for genomeFile in genomes
                           if genomeFile.endswith(('.fna', '.fasta'))]

        genomes = [PrimerEngine.GeneticEntities.Genome(genomeFilePath)
                   for genomeFilePath in genomeFilePaths]

        print("Loaded %i genomes." % len(genomes))

        maxGenomes = 25
        if len(genomes) > maxGenomes:
            print("Discarding genomes, max is %i!" % maxGenomes)

        genomes = genomes[:maxGenomes]
    else:
        genomes = []

    if not genomes:
        print("Fatal: No genomes found!")
        exit(1)

    # APPLY GENOME FEATURES TO BRUTE FORCE MODULE;
    bruteForceSearcher = PrimerEngine.bruteForcePrimerSearch.bruteForceSearcher(genomeFeatures, genomeFilePaths)
    if not bruteForceSearcher.matchedGenome:
        bruteForceSearcher = None

    # -- SETUP OUTPUT DATA STRUCTURES;
    AllLociAmpliconSet = {}
    AllLociPrimerSet = OrderedDict()

    # after this number of tries, we give up on matching primers for the locus.
    RebootLocusTolerance = 13
    matchedPrimerSequences = []

    # ITERATE LOCI;
    for i in range(lociPrimerList.shape[0]):
        locus_info = lociPrimerList.iloc[i]
        locus_name = locus_info["LocusName"]

        # ASSIGN OUTPUT FASTA FILE NAME AND CHECK IF EXISTS;
        outputFastaName = "LOCI_%s.fasta" % locus_name

        outputFastaPath = os.path.join(options.outputPath, outputFastaName)
        if os.path.isfile(outputFastaName):
            print("Skipping locus %s. Already exists..." % locus_name)

        # MAYBE WE WANT TO SKIP GIVEN LOCUS?
        if options.WantedLoci:
            WantedLoci = options.WantedLoci.split(',')
            WantedLoci = [l.strip() for l in WantedLoci]
            if locus_name not in WantedLoci:
                continue

        overallProgress = (i + 1, lociPrimerList.shape[0])

        (LocusAmpliconSet, matchSuccess, primerPair) =\
            PrimerEngine.PrimerDock.matchLocusOnGenomes(
                locus_name,
                locus_info,
                genomes,
                overallProgress,
                bruteForceSearcher=bruteForceSearcher
            )

        if LocusAmpliconSet is not None:
            score = PrimerEngine.ampliconSanity.evaluateSetOfAmplicons(LocusAmpliconSet)
            print("\tAlignment Health = %.2f%%" % score)
            print()
            # record amplicon and primer data;
            writeFastaFile(outputFastaPath, locus_name,
                           LocusAmpliconSet, clonalReference=clonalReference)

            primerPair["AlignmentHealth"] = score
            matchedPrimerSequences.append(primerPair)
            AllLociPrimerSet[locus_name] = matchSuccess
            # print("Bad Amplicon set for %s! Ignoring...." % locus_name)

        else:
            print("WARNING: PrimerDock failure.")

    if matchedPrimerSequences:
        # SHOW AMPLICON DATABASE;
        print(json.dumps(AllLociAmpliconSet, indent=2))

        # BUILD MATCHED PRIMER DATABASE;
        outputFilePath = os.path.join(options.outputPath, "MatchedPrimers.csv")
        data = pd.DataFrame(matchedPrimerSequences,
                            columns=["LocusName",
                                     *PrimerTypes,
                                     "RebootCount",
                                     "AlignmentHealth"])
        data.to_csv(outputFilePath, index=False)

        # Primer Maps on Guide Genome:
        PrimerData = []
        allPrimers = []
        for Locus in AllLociPrimerSet.keys():
            for Primer in AllLociPrimerSet[Locus]:
                row = Primer[0].to_dict(Locus)
                del row["Chromosome"]
                PrimerData.append(row)
                allPrimers.append(Primer)

        outputFilePath = os.path.join(options.outputPath, "PrimerData.csv")
        data = pd.DataFrame(PrimerData)
        data.to_csv(outputFilePath, index=False)

    else:
        print("No regions found, nothing to do.")

    # NOPE
    # MasterGenome = [g for g in genomes if "ME49" in g.name][0]
    # geneGraphs.plotGeneArea(allPrimers, MasterGenome)


if __name__ == "__main__":

    parser = OptionParser()

    parser.add_option("-p", "--plot",
                      dest="PlotArea",
                      action="store_true", default=False)

    parser.add_option("-l",
                      dest="WantedLoci",
                      default="")

    parser.add_option("-i",
                      dest="primerFile")

    parser.add_option("-o",
                      dest="outputPath")

    parser.add_option("-r",
                      "--locusref",
                      dest="LocusReference")

    parser.add_option("-w",
                      "--rewrite",
                      dest="RewriteFasta")

    options, args = parser.parse_args()

    Execute(options)
