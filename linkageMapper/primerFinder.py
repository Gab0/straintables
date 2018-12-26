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
from Bio import SeqIO

import PrimerEngine

from optparse import OptionParser


parser = OptionParser()

parser.add_option("-p", "--plot",
                  dest="PlotArea",
                  action="store_true", default=False)

parser.add_option("-l",
                  dest="WantedLoci",
                  default="")

parser.add_option("-r",
                  dest="primerFile")

parser.add_option("-o",
                  dest="outputPath")

options, args = parser.parse_args()


if __name__ == "__main__":
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
            if PrimerEngine.PrimerDock.evaluateSetOfAmplicons(LocusAmpliconSet):
                # record amplicon and primer data;
                AllLociAmpliconSet[locus_name] = LocusAmpliconSet
                matchedPrimerSequences.append(primerPair)
                AllLociPrimerSet[locus_name] = matchSuccess
            else:
                print("Bad Amplicon set for %s! Ignoring...." % locus_name)

    # SHOW AMPLICON DATABASE;
    print(json.dumps(AllLociAmpliconSet, indent=2))

    # BUILD OUTPUT AMPLICON DATABASE;
    LOCI_SEQUENCE_DATA = []
    data_columns = ["Genome"]

    # TRANPOSE AMPLICON DATABASE!
    TALAS = {}
    for Loci in AllLociAmpliconSet.keys():
        for Genome in AllLociAmpliconSet[Loci].keys():
            if Genome not in TALAS.keys():
                TALAS[Genome] = {}
            TALAS[Genome][Loci] = AllLociAmpliconSet[Loci][Genome]

    AllLociAmpliconSet = TALAS

    for Genome in AllLociAmpliconSet.keys():
        row = {
            "Genome": Genome
        }
        for locus in AllLociAmpliconSet[Genome].keys():
            row[locus] = AllLociAmpliconSet[Genome][locus]
            if locus not in data_columns:
                data_columns.append(locus)
        LOCI_SEQUENCE_DATA.append(row)

    # BUILD OUTPUT AMPLICON DATABASE PATH AND SAVE.
    outputPath = os.path.join(options.outputPath, "Sequences.csv")
    data = pd.DataFrame(LOCI_SEQUENCE_DATA, columns=data_columns)
    data.to_csv(outputPath, index=False)

    # BUILD MATCHED PRIMER DATABASE;
    outputFile = os.path.splitext(os.path.basename(options.primerFile))[0] + "_real.csv"
    outputFilePath = os.path.join(options.outputPath, outputFile)
    data = pd.DataFrame(matchedPrimerSequences,
                        columns=["LocusName", *PrimerTypes, "RebootCount"])
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

