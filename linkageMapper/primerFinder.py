#!/bin/python
import os
import re
import itertools
import json
import random
import numpy as np
import pandas as pd

from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

import geneGraphs
import downloadGene
import primerDockEngine
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



"""

findPrimer:
returns the list of matched primers along with the sequence modification on the primer that resulted in match.

"""


def findPrimer(genome_segment, primer_sequence):
    Primer = Seq(primer_sequence)
    seqVarNames = [
        "Raw Primer",
        "Reverse Complement"
    ]

    sequenceVariations = [
        Primer,
        Primer.reverse_complement(),
    ]

    for n, sequenceVariation in enumerate(sequenceVariations):
        search_seq = str(sequenceVariation).lower()
        Matches = re.finditer(search_seq, genome_segment.sequence)
        Matches = list(Matches)
        if Matches:
            return Matches, seqVarNames[n]

    return [], None


def fetchGeneSequence(geneName, outputFilePath):
    # FETCH PRIMER METHODS. TO BE INTEGRATED;
    DownloadedGene = downloadGene.retrieveGeneSequence(genomeFeatures, geneName)
    if DownloadedGene == None:
        print("ABORT: Gene name not found.")
        exit(1)

    chr_descriptor, location = DownloadedGene
    genomeFilePath = "genomes/TG_ME49.fasta"
    SEQ = downloadGene.fetchSequence(location,
                                     chr_descriptor,
                                     genomeFilePath)

    # Save sequence;
    outputFile = open(outputFilePath, 'w')
    outputFile.write(str(SEQ))
    outputFile.close()


def findPrimerBruteForce(genome, gene_sequence, Reverse=False):
    PRIMER_LENGTH = 20
    SEARCH_STEP = 5

    # FOCUS SEARCH ON A REGION ON THE MIDDLE OF THE GENE SEQUENCE;
    sequenceLength = len(gene_sequence)
    sequenceLengthAim = 1500
    if sequenceLength > sequenceLengthAim:
        sequenceLengthBounds = (
            sequenceLength // 2 - sequenceLengthAim // 2,
            sequenceLength // 2 + sequenceLengthAim // 2
        )
        gene_sequence = gene_sequence[sequenceLengthBounds[0]:sequenceLengthBounds[1]]

    if Reverse:
        Indexes = range(len(gene_sequence) - PRIMER_LENGTH, 0, -SEARCH_STEP)
    else:
        Indexes = range(0, len(gene_sequence), SEARCH_STEP)

    foundPrimers = []
    for s in Indexes:
        primer_sequence = gene_sequence[s:s + PRIMER_LENGTH]
        for c, _chr in enumerate(genome):
            matches, sequenceVariationName = findPrimer(_chr, primer_sequence)
            if len(matches) > 1:
                print("Leak.")
                continue
            if matches:
                print(matches[0][0].upper())
                print(sequenceVariationName)
                foundPrimers.append(matches[0])
                if len(foundPrimers) > BFPCL:
                    return foundPrimers

        if (s <= (len(gene_sequence) // 5)) == Reverse:
            break

    return foundPrimers


def matchPrimerOnGenome(genome, PrimerSequence, PrimerType):
    matchedPrimers = []

    for chromosome in genome:
        primerMatches, sequenceVariationName = findPrimer(
            chromosome,
            PrimerSequence
        )
        if primerMatches:
            print("\t@%s" % sequenceVariationName)

            # TBD: what to do about this?
            if len(primerMatches) > 1:
                print("Primer overflow!")

            match = primerMatches[0]
            match_position = match.start()

            # SEARCH REGION ON ANNOTATED CHROMOSOMES;
            print("Found at chromosome %s" % chromosome.name)
            print("pos %i of %i" % (match.start(), chromosome.length))

            # RECORD PRIMER;
            matchedPrimer = primerDockEngine.primerMatch(match,
                                                         PrimerType,
                                                         chromosome,
                                                         PrimerSequence
            )

            matchedPrimers.append(matchedPrimer)
            if len(matchedPrimers) > BFPCL:
                return matchedPrimers

    return matchedPrimers


def launchBruteForcePrimerSearch(locus_name, chromosomes, Reverse):

    # BRUTE FORCE PRIMER FINDER OPERATIONS;
    geneSequenceFile = "%s.fasta" % locus_name
    geneSequenceFilePath = os.path.join("PrimerSources", geneSequenceFile)

    if not os.path.isfile(geneSequenceFilePath):
        # Fetch gene sequence;
        fetchGeneSequence(locus_name, geneSequenceFilePath)

    geneSequenceRaw = open(geneSequenceFilePath).read()

    geneSequence = geneSequenceRaw.split("\n")
    if ">" in geneSequence[0]:
        geneSequence = geneSequence[1:]
    geneSequence = "".join(geneSequence).lower()

    foundPrimers = findPrimerBruteForce(chromosomes, geneSequence, Reverse=Reverse)
    if foundPrimers:
        print("Brute force forward primer -> %s" % foundPrimers)

    resultingPrimers = [p[0].upper() for p in foundPrimers] 

    return resultingPrimers


def searchPrimerPairOnGenome(locusName, primerPair, genome):

    # INIT SUPPORT LIST;
    matchedPrimers = {}

    matchSuccess = [True, True]
    # ITERATE THRU PRIMER TYPES;
    for PT, PrimerType in enumerate(PrimerTypes):

        # screen info purposes;
        print(PrimerType)

        # fetch primer sequence;
        queryPrimer = primerPair[PrimerType]
        print(queryPrimer)
        matchedPrimers[PrimerType] = matchPrimerOnGenome(genome,
                                                         queryPrimer,
                                                         PrimerType)

        if not matchedPrimers[PrimerType]:
            print("No match...")
            matchSuccess[PT] = False

    for PT, PrimerType in enumerate(PrimerTypes):
        # CHECK FOR PRIMER LEAK;
        PrimerType = PrimerTypes[PT]

        if len(matchedPrimers[PrimerType]) > 1:
            print("Primer leak... trying to fix.")

            # TRY TO FIX PRIMER LEAK;
            # outright delete matches that doesn't have a pair in the same chromosome;
            # if PrimerTypes[1 - PT] in matchedPrimers.keys():
            opponentIndexes = [P.chr_index for P in matchedPrimers[PrimerTypes[1 - PT]]]
            for p, P in enumerate(matchedPrimers[PrimerType]):
                if P.chr_index not in opponentIndexes:
                    matchedPrimers[PrimerType][p] = None
            matchedPrimers[PrimerType] = [p for p in matchedPrimers[PrimerType] if p]

            # might be unnecessary
            matchedPrimers[PrimerType] = sorted(matchedPrimers[PrimerType],
                                                key=lambda x: x.chr_length,
                                                reverse=True)[:1]

            # PRIMER LEAK IS UNAVOIDABLE;
            if len(matchedPrimers[PrimerType]) > 1:
                matchSuccess[PT] = False
        print()

    print("\n\n")

    # RETRIEVE INTERPRIMER SEQUENCES;
    matchCount = sum([len(matchedPrimers[PrimerType])
                      for PrimerType in PrimerTypes])

    # TWO MATCHES... IDEAL SCENARIO;
    if matchCount == 2:
        Primers = [matchedPrimers[PrimerTypes[i]][0] for i in range(2)]
        # if anything goes wrong while building the amplicon, the match fails.
        # Amplicon may rase errors deliberately when match result is not optimal.
        try:
            amplicon = primerDockEngine.Amplicon(genome, Primers[0], Primers[1])
        except ValueError:
            return "", [False, False]

        return amplicon.Sequence, (matchedPrimers[PrimerTypes[0]], matchedPrimers[PrimerTypes[1]])
    else:
        return "", matchSuccess


def validatePrimer(V):
    if type(V) == np.float64:
        return False
    if type(V) == float:
        return False
    elif not V:
        return False
    else:
        return True


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
    genomes = [primerDockEngine.Genome(genomeFilePath)
               for genomeFilePath in genomeFilePaths]

    AllLociAmpliconSet = {}
    AllLociPrimerSet = OrderedDict()
    AllLociSuccessEvaluation = OrderedDict()

    # brute force primer count limiter;
    BFPCL = 36

    # after this number of tries, we give up on matching primers for the locus.
    RebootLocusTolerance = 150
    matchedPrimerSequences = []
    # ITERATE LOCI;
    
    for i in range(lociPrimerList.shape[0]):
        locus_info = lociPrimerList.iloc[i]
        locus_name = locus_info["LocusName"]

        LocusAmpliconSet = {}

        if options.WantedLoci:
            WantedLoci = options.WantedLoci.split(',')
            WantedLoci = [l.strip() for l in WantedLoci]
            if locus_name not in WantedLoci:
                continue

        def initPrimerHolder():
            return {k: [] for k in locus_info.keys()}

        primerTrash = initPrimerHolder()
        testablePrimers = initPrimerHolder()


        # load primer pair data from user-made Primer file;
        primerPair = dict(locus_info)

        RebootCount = 0

        # ITERATE GENOMES UNTIL SUCCESS;
        for Genome in itertools.cycle(genomes):
            # print("Genome: %s --\n" % genome)

            HEADER_INFO = (i + 1, lociPrimerList.shape[0],
                           RebootCount + 1, locus_name, Genome.name)
            print(">>> Locus %i of %i | run number %i ->  Searching %s on %s" % HEADER_INFO)
            # primer sequence may baae uknown;
            matchSuccess = [True, True]
            for PT, PrimerType in enumerate(PrimerTypes):
                if not validatePrimer(primerPair[PrimerType]):
                    matchSuccess[PT] = False

            AmpliconSequence = ""
            print("\n")

            if all(matchSuccess):
                print("Searching sequence for locus %s" % locus_name)

                AmpliconSequence, matchSuccess =\
                    searchPrimerPairOnGenome(locus_name, primerPair, Genome)

                print(matchSuccess)

            # FAILURE ON MATCHING A PRIMER?
            if not all(matchSuccess):
                for PT, match in enumerate(matchSuccess):
                    PrimerType = PrimerTypes[PT]
                    if not match:
                        print("Resetting %s!" % PrimerType)

                        # delete current primer;
                        if validatePrimer(primerPair[PrimerType]):
                            primerTrash[PrimerType].append(primerPair[PrimerType])

                        primerPair[PrimerType] = None

                        # RESET ALL MATCHES!
                        LocusAmpliconSet = {}
                        LocusPrimerSet = {}

                        RebootCount += 1
                        # brute force for the problematic genome!
                        print("Searching new primer on gene sequence...")

                        # MANAGE PRIMERS ON HOLSTER;
                        if not testablePrimers[PrimerType]:
                            testablePrimers[PrimerType] =\
                                launchBruteForcePrimerSearch(locus_name, Genome, PT)

                        while testablePrimers[PrimerType]:
                            print("Fetching new primer from reserve...")
                            random.shuffle(testablePrimers[PrimerType])
                            primerPair[PrimerType] = testablePrimers[PrimerType].pop()
                            if primerPair[PrimerType] not in primerTrash[PrimerType]:
                                break
            else:
                # SUCCESS.
                if len(AmpliconSequence) > 100:
                    print("Found amplicon of length %i" % len(AmpliconSequence))
                    LocusAmpliconSet[Genome.name] = AmpliconSequence
                    # AMPLICON IS TOO SHORT;
                else:
                    print("Resetting all primers, amplicon too short.")
                    for s in range(2):
                        if validatePrimer(primerPair[PrimerTypes[s]]):
                            primerTrash[PrimerTypes[s]].append(primerPair[PrimerTypes[s]])
                        primerPair[PrimerTypes[s]] = None

            # CHECK LOCUS PROGRESS ACROSS ALL GENOMES;
            progressMark = len(list(LocusAmpliconSet.keys())) / len(genomes)
            print("%s match progress: %.2f%%" % (locus_name, progressMark * 100))
            print()

            # IS LOCUS DONE?
            if progressMark == 1:
                print(">>> Matched %s on all genomes." % locus_name)
                print()
                print()

                # export amplicon and primer data;
                AllLociAmpliconSet[locus_name] = LocusAmpliconSet
                matchedPrimerSequences.append(primerPair)
                AllLociPrimerSet[locus_name] = matchSuccess

                AllLociSuccessEvaluation[locus_name] = RebootCount
                # exit loop...
                break

            # IF IT FAILS ENOUGH, SKIP LOCUS;
            elif RebootCount > RebootLocusTolerance:
                break

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
    data = pd.DataFrame(matchedPrimerSequences, columns=["LocusName"] + PrimerTypes)
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

    print(json.dumps(AllLociSuccessEvaluation, indent=2))
