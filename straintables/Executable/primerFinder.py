#!/bin/python
import os
import numpy as np
import pandas as pd

from collections import OrderedDict

from Bio import Seq, SeqIO

import argparse

from straintables import PrimerEngine, InputFile, OutputFile, Definitions
from straintables.Database import annotationManager, genomeManager


def writeFastaFile(outputPath,
                   locusName,
                   locusSequences):
    fastaSequences = []

    for genome in locusSequences.keys():
        sequence = SeqIO.SeqRecord(Seq.Seq(locusSequences[genome]),
                                   id=genome,
                                   name=genome,
                                   description="")
        fastaSequences.append(sequence)

    with open(outputPath, "w") as output_handle:
        SeqIO.write(fastaSequences, output_handle, "fasta")


def Execute(options):

    # CHECK DECLARATION OF PRIMER FILE;
    if not options.PrimerFile:
        print("FATAL: No primer file specified.")
        exit(1)

    print("\nAnnotation feature type for automatic primer search is |%s|." %
          options.wantedFeatureType)
    print("\t (allowed options: CDS, gene, mRNA)")

    # -- CHECK ANNOTATED FEATURES; -- useles step maybe
    featureFolderPath =\
        os.path.join(options.SourceDataDirectory, "annotations")

    if os.path.isdir(featureFolderPath):
        genomeFeatureFiles = [
            os.path.join(featureFolderPath, File)
            for File in os.listdir(featureFolderPath)
            if not File.startswith(".")
        ]

    else:
        genomeFeatureFiles = []

    # CHECK GENOME FEATURES FILE EXISTENCE;
    if not genomeFeatureFiles:
        print("Fatal: No genbank features file found.")
        exit(1)

    # -- LOAD USER DEFINED PRIMERS;
    lociPrimerList = InputFile.loadPrimerList(options.PrimerFile)

    # LOAD GENOMES;
    genomeFilePaths = genomeManager.readGenomeFolder(
        os.path.join(options.SourceDataDirectory, "genomes")
    )
    genomes = [PrimerEngine.GeneticEntities.Genome(genomeFilePath)
               for genomeFilePath in genomeFilePaths]

    print("Loaded %i genomes." % len(genomes))

    maxGenomes = 100
    if len(genomes) > maxGenomes:
        print("Discarding genomes, max is %i!" % maxGenomes)

        genomes = genomes[:maxGenomes]

    if not genomes:
        print("Fatal: No genomes found!")
        exit(1)

    if len(genomes) < 4:
        print("Fatal: need at least 4 genomes to proper execute the analysis,")
        print("\tgot only %i." % len(genomes))
        exit(1)

    # Avoid repeating genome names;
    genomeNames = []
    for genome in genomes:
        while genome.name in genomeNames:
            genome.name += "+"
        genomeNames.append(genome.name)

    # Initialize brute force prime searcher on top of chosen annotation file;
    annotationFilePath, genomeFeatures =\
        annotationManager.loadAnnotation(featureFolderPath)

    if annotationFilePath:
        bruteForceSearcher =\
            PrimerEngine.PrimerDesign.BruteForcePrimerSearcher(
                genomeFeatures,
                genomeFilePaths,
                wantedFeatureType=options.wantedFeatureType,
                FindPCRViablePrimers=options.RealPrimers,
                PrimerLength=options.PrimerLength,
                AmpliconMinimumLength=options.MinAmpliconLength,
                AmpliconMaximumLength=options.MaxAmpliconLength,
                PrimerAllowedUncertainty=options.PrimerUncertainty
            )

    if not annotationFilePath:
        # or not bruteForceSearcher.matchedGenome:(matchedGenome is deprecated)
        bruteForceSearcher = None

    # -- SETUP OUTPUT DATA STRUCTURES;
    AllLociPrimerSet = OrderedDict()

    matchedPrimerSequences = []

    print("\n")

    GenomeFailureReport =\
        OutputFile.DockFailureReport(options.WorkingDirectory)

    # ITERATE LOCI: Main Loop;
    for i in range(lociPrimerList.shape[0]):
        locus_info = lociPrimerList.iloc[i]
        locus_name = locus_info["LocusName"]

        # ASSIGN OUTPUT FASTA FILE NAME AND CHECK IF EXISTS;
        outputFastaName = "%s%s.fasta" % (
            Definitions.FastaRegionPrefix, locus_name)

        outputFastaPath = os.path.join(
            options.WorkingDirectory, outputFastaName)

        print("Fasta file: %s" % outputFastaPath)

        if os.path.isfile(outputFastaPath):
            print("Skipping locus %s. Already exists..." % locus_name)
            continue

        # MAYBE WE WANT TO SKIP GIVEN LOCUS?
        if options.WantedLoci:
            WantedLoci = options.WantedLoci.split(',')
            WantedLoci = [l.strip() for l in WantedLoci]
            if locus_name not in WantedLoci:
                continue

        overallProgress = (i + 1, lociPrimerList.shape[0])

        RegionMatchResult =\
            PrimerEngine.PrimerDock.matchLocusOnGenomes(
                locus_name,
                locus_info,
                genomes,
                overallProgress,
                rebootTolerance=options.rebootTolerance,
                allowN=options.AllowUnknownBaseInAmplicon,
                bruteForceSearcher=bruteForceSearcher
            )

        FailureType = PrimerEngine.PrimerDock.RegionMatchFailure
        if type(RegionMatchResult) == FailureType:
            GenomeFailureReport.content[locus_name] =\
                RegionMatchResult.FailedGenomes
            continue

        # -- Additional region statistics;
        if RegionMatchResult.LocusAmpliconSet is not None:
            # AlignmentHealth.
            score = PrimerEngine.ampliconSanity.evaluateSetOfAmplicons(
                RegionMatchResult.LocusAmpliconSet)

            print("\tAlignment Health = %.2f%%" % score)
            print()
            # record amplicon and primer data;
            writeFastaFile(outputFastaPath, locus_name,
                           RegionMatchResult.LocusAmpliconSet)

            primerPair = {
                P.label: P.sequence
                for P in RegionMatchResult.MatchedPrimers
            }

            primerPair["LocusName"] = locus_name

            primerPair["AlignmentHealth"] = score

            RegionLengths = [
                len(r)
                for r in RegionMatchResult.LocusAmpliconSet
            ]

            primerPair["MeanLength"] = np.mean(RegionLengths)
            primerPair["StdLength"] = np.std(RegionLengths)
            primerPair["Chromosome"] = RegionMatchResult.chr_identifier
            primerPair["StartPosition"] =\
                RegionMatchResult.MatchedPrimers[0].position.start()

            # Append region data;
            matchedPrimerSequences.append(primerPair)
            AllLociPrimerSet[locus_name] = RegionMatchResult.MatchedPrimers
            # print("Bad Amplicon set for %s! Ignoring...." % locus_name)
        else:
            print("WARNING: PrimerDock failure.")

    if matchedPrimerSequences:
        # SHOW AMPLICON DATABASE;

        # BUILD MATCHED PRIMER DATABASE;
        MatchedRegions = OutputFile.MatchedRegions(options.WorkingDirectory)
        MatchedRegions.add(matchedPrimerSequences)
        MatchedRegions.write()

        # Primer Maps on Guide Genome:
        PrimerData = []
        allPrimers = []

        for Locus in AllLociPrimerSet.keys():
            for Primer in AllLociPrimerSet[Locus]:
                row = Primer.to_dict(Locus)

                del row["Chromosome"]
                PrimerData.append(row)
                allPrimers.append(Primer)

        # -- SAVE PRIMER DATA FILE;
        fPrimerData = OutputFile.PrimerData(options.WorkingDirectory)
        fPrimerData.content = pd.DataFrame(PrimerData)
        fPrimerData.write()

        # BUILD INFORMATION FILE;
        information = OutputFile.AnalysisInformation(options.WorkingDirectory)
        information.content = {
            "annotation": os.path.realpath(annotationFilePath),
            }
        information.content.update(options.__dict__)
        information.write()

        # -- WRITE GENOME FAILURE REPORT IF ANY REGION FAILED;
        if GenomeFailureReport.content:
            GenomeFailureReport.write()
    else:
        print("No regions found, nothing to do.")

    # NOPE
    # MasterGenome = [g for g in genomes if "ME49" in g.name][0]
    # geneGraphs.plotGeneArea(allPrimers, MasterGenome)
    return matchedPrimerSequences


def parse_arguments(parser):

    parser.add_argument("--plot",
                        dest="PlotArea",
                        action="store_true", default=False)

    parser.add_argument("-l",
                        dest="WantedLoci",
                        default="")

    parser.add_argument("-p",
                        dest="PrimerFile")

    parser.add_argument("-o",
                        dest="outputPath")

    parser.add_argument("-r",
                        "--locusref",
                        dest="LocusReference")

    parser.add_argument("-w",
                        "--rewrite",
                        dest="RewriteFasta")

    parser.add_argument("-t",
                        "--feature",
                        dest="wantedFeatureType",
                        default="gene",
                        help="Wanted feature type at annotation " +
                        "for primer design. [gene|CDS|mRNA]"),

    parser.add_argument("--pcr",
                        dest="RealPrimers",
                        action="store_true",
                        help="Search only for primers " +
                        "that are viable in real PCR reactions.")

    parser.add_argument("-b", "--reboot",
                        dest="rebootTolerance",
                        type=int,
                        default=20)

    parser.add_argument("--maxalignlen",
                        dest="MaxAmpliconLength",
                        type=int,
                        default=1200,
                        help="Maximum sequence length for extracted regions.")

    parser.add_argument("--minalignlen",
                        dest="MinAmpliconLength",
                        type=int,
                        default=400,
                        help="Minimum sequence length for extracted regions.")

    parser.add_argument("--allow-n",
                        dest="AllowUnknownBaseInAmplicon",
                        action="store_true")

    parser.add_argument("-d", "--dir", dest="WorkingDirectory")

    parser.add_argument("--datadir", dest="SourceDataDirectory", default="")

    parser.add_argument("--primerlength", dest="PrimerLength",
                        default=20, type=int)

    parser.add_argument("--primer-uncertainty", dest="PrimerUncertainty",
                        default=0.0, type=float,
                        help="Percent of primer letters that can be set to match any base, by using the 'letter' '.'.")
    return parser


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    options = parse_arguments(parser).parse_args()
    Execute(options)
