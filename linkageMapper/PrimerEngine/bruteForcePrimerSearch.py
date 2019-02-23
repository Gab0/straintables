#!/bin/python

import os

from . import PrimerDock
from Bio import SeqIO

"""

This finds a feature compatible to the query gene name, in a chromosome feature table.

returns: (Name of chromosome, position of sequence inside chromosome)
"""


def retrieveGeneSequence(genomeFeatures, geneName):
    for g, FeatureGroup in enumerate(genomeFeatures):
        for feature in FeatureGroup.features:
            if feature.type == "gene":
                MATCH = False
                if "gene" in feature.qualifiers.keys():
                    if geneName in feature.qualifiers['gene']:
                        MATCH = True
                else:
                    if geneName in feature.qualifiers['locus_tag']:
                        MATCH = True
                if MATCH:
                    return FeatureGroup.description, feature.location


def locateAndFetchSequence(location, chr_descriptor, genomeFilePath):
    genome = list(SeqIO.parse(genomeFilePath, format="fasta"))

    wantedDescriptors = [chr_descriptor, "complete genome"]
    for c, Chromosome in enumerate(genome):
        for Descriptor in wantedDescriptors:
            if Descriptor in Chromosome.description:
                Sequence = Chromosome.seq[location.start.position:location.end.position]
                if location.strand == -1:
                    Sequence = Sequence.reverse_complement()
                return Sequence
        # print(dir(Chromosome))
        # print(Chromosome.name)
        # print(Chromosome.id)
        # print(Chromosome.description)


def fetchGeneSequence(geneName, outputFilePath):

    # FETCH PRIMER METHODS. TO BE INTEGRATED;
    geneSearchResult = retrieveGeneSequence(genomeFeatures, geneName)

    if geneSearchResult is None:
        print("Aborting brute force primer search: Gene name not found.")
        return

    chr_descriptor, location = geneSearchResult

    AnnotationDescription = genomeFeatures[0].description.split(" ")

    # print(AnnotationDescription)

    matchingAnnotationGenome = None
    for k in genomeFilePaths:
        for m in AnnotationDescription:
            if m.lower() in k:
                matchingAnnotationGenome = k

    if matchingAnnotationGenome is None:
        print("No genome matching annotation!")
        exit(1)

    SEQ = locateAndFetchSequence(location,
                                 chr_descriptor,
                                 matchingAnnotationGenome)

    if not SEQ:
        print("Error: Failure on feching brute force sequence.")
        print("genomePath: %s" % matchingAnnotationGenome)
        print("chromosome descripor: %s" % chr_descriptor)
        print("location: %s" % location)
        exit(1)

    # Save sequence;
    outputFile = open(outputFilePath, 'w')
    outputFile.write(str(SEQ))
    outputFile.close()


def launchBruteForcePrimerSearch(locus_name, chromosomes, Reverse):

    # BRUTE FORCE PRIMER FINDER OPERATIONS;
    geneSequenceFile = "%s.fasta" % locus_name
    geneSequenceFilePath = os.path.join("PrimerSources", geneSequenceFile)

    if not os.path.isfile(geneSequenceFilePath):
        # Fetch gene sequence;
        fetchGeneSequence(locus_name, geneSequenceFilePath)

    if os.path.isfile(geneSequenceFilePath):
        geneSequenceRaw = open(geneSequenceFilePath).read()
    else:
        return None

    geneSequence = geneSequenceRaw.split("\n")
    if ">" in geneSequence[0]:
        geneSequence = geneSequence[1:]
    geneSequence = "".join(geneSequence).lower()

    foundPrimers =\
        findPrimerBruteForce(chromosomes,
                             geneSequence,
                             Reverse=Reverse
        )

    if foundPrimers:
        print("Brute force forward primer search returns %i primers." % len(foundPrimers))

    resultingPrimers = [p[0].upper() for p in foundPrimers] 

    return resultingPrimers


def findPrimerBruteForce(genome, gene_sequence,
                         Reverse=False, maximumPrimerCount=36, Verbose=False):
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
            matches, sequenceVariationName = PrimerDock.findPrimer(_chr, primer_sequence)
            if len(matches) > 1:
                if Verbose:
                    print("Leak.")
                continue
            if matches:
                if Verbose:
                    print(matches[0][0].upper())
                    print(sequenceVariationName)
                foundPrimers.append(matches[0])
                if len(foundPrimers) > maximumPrimerCount:
                    return foundPrimers

        if (s <= (len(gene_sequence) // 5)) == Reverse:
            break

    return foundPrimers


