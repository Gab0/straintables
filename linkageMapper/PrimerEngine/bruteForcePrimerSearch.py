#!/bin/python

import os

from . import PrimerDock
from Bio import SeqIO


def retrieveGeneSequence(genomeFeatures, geneName):
    for g, FeatureGroup in enumerate(genomeFeatures):
        for feature in FeatureGroup.features:
            A = FeatureGroup
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


def fetchSequence(location, chr_descriptor, genomeFilePath):
    genome = list(SeqIO.parse(genomeFilePath, format="fasta"))

    for c, Chromosome in enumerate(genome):
        if chr_descriptor in Chromosome.description:
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
    DownloadedGene = retrieveGeneSequence(genomeFeatures, geneName)
    if DownloadedGene is None:
        print("ABORT: Gene name not found.")
        exit(1)

    chr_descriptor, location = DownloadedGene
    genomeFilePath = "genomes/TG_ME49.fasta"
    SEQ = fetchSequence(location,
                        chr_descriptor,
                        genomeFilePath)

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

    geneSequenceRaw = open(geneSequenceFilePath).read()

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
                         Reverse=False, maximumPrimerCount=36):
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
                print("Leak.")
                continue
            if matches:
                print(matches[0][0].upper())
                print(sequenceVariationName)
                foundPrimers.append(matches[0])
                if len(foundPrimers) > maximumPrimerCount:
                    return foundPrimers

        if (s <= (len(gene_sequence) // 5)) == Reverse:
            break

    return foundPrimers


