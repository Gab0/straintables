#!/bin/python

import os

from . import PrimerDock
from Bio import SeqIO

"""

This finds a feature compatible to the query gene name, in a chromosome feature table.

returns: (Name of chromosome, position of sequence inside chromosome)
"""


class bruteForceSearcher():
    def __init__(self, genomeFeatures, genomeFilePaths):
        self.genomeFeatures = genomeFeatures
        self.matchedGenome = self.locateMatchingGenome(genomeFilePaths)
        if self.matchedGenome is None:
            print()
            print("Warning: automatic primer search disabled. No matching genome found!")
            print()
            return None

    def locateMatchingGenome(self, genomeFilePaths, Verbose=True):

        AnnotationDescription = self.genomeFeatures[0].description
        if Verbose:
            print("Searching a genome that matches the annotation:")
            print(AnnotationDescription)

        matchingGenomeFilePath = None
        for k in genomeFilePaths:
            strippedGenomeName = os.path.splitext(os.path.split(k)[-1])[0]
            if strippedGenomeName in AnnotationDescription:
                matchingGenomeFilePath = k
                print(strippedGenomeName)

        if matchingGenomeFilePath is None:
            print("No genome matching annotation!")
            return None
        else:
            print("Found matching for automatic primer search: %s" % matchingGenomeFilePath)

        genome = list(SeqIO.parse(matchingGenomeFilePath, format="fasta"))

        return genome

    def retrieveGeneSequence(self, geneName):
        for g, FeatureGroup in enumerate(self.genomeFeatures):
            #print(FeatureGroup)
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
        print("Warning: Gene %s not found." % geneName)
        #exit(1)
    def locateAndFetchSequence(self, location, chr_descriptor):

        wantedDescriptors = [chr_descriptor, "complete genome"]
        for c, Chromosome in enumerate(self.matchedGenome):
            for Descriptor in wantedDescriptors:
                print("Fetching primers from %s..." % Descriptor)
                if Descriptor in Chromosome.description:
                    Sequence = Chromosome.seq[location.start.position:location.end.position]
                    if location.strand == -1:
                        Sequence = Sequence.reverse_complement()
                    return Sequence

    def fetchGeneSequence(self, geneName, outputFilePath):

        # FETCH PRIMER METHODS. TO BE INTEGRATED;
        geneSearchResult = self.retrieveGeneSequence(geneName)

        if geneSearchResult is None:
            print("Aborting brute force primer search: Gene name not found.")
            return

        chr_descriptor, location = geneSearchResult


        # print(AnnotationDescription)


        SEQ = self.locateAndFetchSequence(location,
                                          chr_descriptor
                                          )

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


    def launchBruteForcePrimerSearch(self, locus_name, chromosomes, Reverse):

        # BRUTE FORCE PRIMER FINDER OPERATIONS;
        geneSequenceFile = "%s.fasta" % locus_name

        PrimerSourcesDirectory = "PrimerSources"
        if not os.path.isdir(PrimerSourcesDirectory):
            os.mkdir(PrimerSourcesDirectory)
            
        geneSequenceFilePath = os.path.join(PrimerSourcesDirectory, geneSequenceFile)

        if not os.path.isfile(geneSequenceFilePath):
            # Fetch gene sequence;
            self.fetchGeneSequence(locus_name, geneSequenceFilePath)

        if os.path.isfile(geneSequenceFilePath):
            geneSequenceRaw = open(geneSequenceFilePath).read()
        else:
            return None

        geneSequence = geneSequenceRaw.split("\n")
        if ">" in geneSequence[0]:
            geneSequence = geneSequence[1:]
        geneSequence = "".join(geneSequence).lower()

        foundPrimers =\
            self.findPrimerBruteForce(chromosomes,
                                      geneSequence,
                                      Reverse=Reverse
            )

        if foundPrimers:
            print("Brute force forward primer search returns %i primers." % len(foundPrimers))

        resultingPrimers = [p[0].upper() for p in foundPrimers] 

        return resultingPrimers

    def findPrimerBruteForce(self,
                             genome,
                             gene_sequence,
                             Reverse=False,
                             maximumPrimerCount=36,
                             Verbose=False):
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


