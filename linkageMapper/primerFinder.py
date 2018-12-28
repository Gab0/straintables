#!/bin/python
import os
import re
import itertools
import json
import random
import numpy as np
import pandas as pd

from collections import OrderedDict

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import Seq, SeqIO

import PrimerEngine

from optparse import OptionParser


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


def writeFastaFile(outputPath, locusName, locusSequences):
    fastaSequences = []

    for genome in locusSequences.keys():
        Name = genome + ".fasta"
        try:
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


if __name__ == "__main__":

    # LOAD CLONAL TYPE LOCUS INFORMATION (Su et al.);
    clonalReference = clonalTypeReference()

    # CHECK DECLARATION OF PRIMER FILE;
    if not options.primerFile:
        print("FATAL: No primer file specified.")
        exit(1)

    PrimerTypes = ["ForwardPrimer", "ReversePrimer"]

    # LOAD GENOME FEATURES;
    featureFolderPath = "annotations"
    genomeFeatureFiles = [
        os.path.join(featureFolderPath, File)
        for File in os.listdir(featureFolderPath)
        if not File.startswith(".")
    ]

    genomeFeatures = [
        SeqIO.read(File, "genbank")
        for File in genomeFeatureFiles
    ]

    PrimerEngine.bruteForcePrimerSearch.genomeFeatures = genomeFeatures

    genomeFeaturesChromosomes =\
        [chromosome.features[0].qualifiers['chromosome'][0].upper()
         for chromosome in genomeFeatures]

    print(genomeFeaturesChromosomes)

    # load source data;
    lociPrimerList = pd.read_csv(options.primerFile)

    # LOAD GENOMES;
    genomeDirectory = "genomes"
    genomes = os.listdir(genomeDirectory)
    genomeFilePaths = [os.path.join(genomeDirectory, genomeFile)
                       for genomeFile in genomes
                       if genomeFile.endswith('.fasta')]
    genomes = [PrimerEngine.GeneticEntities.Genome(genomeFilePath)
               for genomeFilePath in genomeFilePaths]

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
        LocusAmpliconSet, matchSuccess, primerPair =\
            PrimerEngine.PrimerDock.matchLocusOnGenomes(
                locus_name,
                locus_info,
                genomes,
                overallProgress
            )

        if LocusAmpliconSet is not None:
            score = PrimerEngine.PrimerDock.evaluateSetOfAmplicons(LocusAmpliconSet)
            # record amplicon and primer data;
            # AllLociAmpliconSet[locus_name] = LocusAmpliconSet
            writeFastaFile(outputFastaPath, locus_name, LocusAmpliconSet)

            primerPair["AlignmentHealth"] = score
            matchedPrimerSequences.append(primerPair)
            AllLociPrimerSet[locus_name] = matchSuccess
            # print("Bad Amplicon set for %s! Ignoring...." % locus_name)

    # SHOW AMPLICON DATABASE;
    print(json.dumps(AllLociAmpliconSet, indent=2))

    # BUILD MATCHED PRIMER DATABASE;
    outputFilePath = os.path.join(options.outputPath, "MatchedPrimers.csv")
    data = pd.DataFrame(matchedPrimerSequences,
                        columns=["LocusName",
                                 *PrimerTypes,
                                 "RebootCount",
                                 "Quality"])
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

    # NOPE
    MasterGenome = [g for g in genomes if "ME49" in g.name][0]
    # geneGraphs.plotGeneArea(allPrimers, MasterGenome)

