#!/bin/python

import os
import random

from . import PrimerDock, RealPrimers
from Bio import SeqIO, Seq

from straintables.Database.StrainNames import fetchStrainName


"""

This finds a feature compatible to the query gene name,
in a chromosome feature table.

returns: (Name of chromosome, position of sequence inside chromosome)
"""


class BruteForcePrimerSearcher():
    def __init__(self,
                 genomeFeatures,
                 genomeFilePaths,
                 wantedFeatureType="CDS",
                 PrimerLength=20,
                 FindPCRViablePrimers=False,
                 AmpliconMinimumLength=400,
                 AmpliconMaximumLength=1200,
                 PrimerAllowedUncertainty=0.0):

        assert(wantedFeatureType in ["gene", "mRNA", "CDS"])

        self.genomeFeatures = genomeFeatures

        # DEPRECATED
        self.matchedGenome = None
        # self.matchedGenome = self.locateMatchingGenome(genomeFilePaths)

        self.wantedFeatureType = wantedFeatureType
        self.PrimerLength = PrimerLength
        self.FindPCRViablePrimers = FindPCRViablePrimers

        self.AmpliconMaximumLength = AmpliconMaximumLength
        self.AmpliconMinimumLength = AmpliconMinimumLength

        self.PrimerAllowedUncertainty = PrimerAllowedUncertainty
        # if self.matchedGenome is None:
        #    print()
        #    print("Warning: automatic primer search disabled.")
        #    print("\tNo matching genome found!")
        #    print()
        #    return None

    def locateMatchingGenome(self, genomeFilePaths, Verbose=False):
        AnnotationDescriptor = self.genomeFeatures[0].description

        matchingGenomeFilePath = None
        annotationStrain = fetchStrainName(AnnotationDescriptor)

        print("\nSearching a genome that matches the annotation..."
              "(strain: %s)" % annotationStrain)

        if Verbose:
            print(AnnotationDescriptor)

        # -- SEARCH BY ANNOTATION INFORMATION;
        for genomePath in genomeFilePaths:
            features = list(SeqIO.parse(genomePath, format="fasta"))
            GenomeDescriptor = features[0].description
            if Verbose:
                print(">%s" % GenomeDescriptor)

            strain = fetchStrainName(GenomeDescriptor)
            if strain and strain == annotationStrain:
                matchingGenomeFilePath = genomePath
                matchingGenomeDescriptor = GenomeDescriptor
                matchingStrain = strain

        if matchingGenomeFilePath is None:
            print("No genome matching annotation!")
            return None
        else:
            print("Found matching genome to annotation, "
                  "for automatic primer search: %s" % matchingGenomeFilePath)
            print("Matching genome descriptor: %s" % matchingGenomeDescriptor)
            print("Detected genome strain: %s" % matchingStrain)

        genome = list(SeqIO.parse(matchingGenomeFilePath, format="fasta"))
        return genome

    def retrieveGeneLocation(self, geneName, wantedFeatureType="CDS"):

        for g, FeatureGroup in enumerate(self.genomeFeatures):
            for feature in FeatureGroup.features:
                if feature.type == wantedFeatureType:
                    MATCH = False
                    if "gene" in feature.qualifiers.keys():
                        if geneName in feature.qualifiers['gene']:
                            MATCH = True
                    else:
                        if geneName in feature.qualifiers['locus_tag']:
                            MATCH = True
                    if MATCH:
                        return FeatureGroup, feature.location

        print("Warning: Gene %s not found." % geneName)

    def locateAndFetchSequence(self, Sequence, position):
        seq = Sequence[position.start.position:position.end.position]
        if position.strand == -1:
            seq = seq.reverse_complement()
        return seq

    def locateAndFetchSequenceOnGenomeFile(self, position, chr_descriptor):
        wantedDescriptors = [chr_descriptor, "complete genome"]
        if not self.matchedGenome:
            print("No matching genome to find gene sequence.")
            return ""
        for c, Chromosome in enumerate(self.matchedGenome):
            for Descriptor in wantedDescriptors:
                print("Fetching primers from %s..." % Descriptor)
                if Descriptor in Chromosome.description:
                    Sequence = Chromosome.seq[
                        position.start.position:position.end.position
                    ]

                    if position.strand == -1:
                        Sequence = Sequence.reverse_complement()
                    return Sequence

    def fetchGeneSequence(self, geneName):

        # -- FETCH PRIMER METHODS.
        geneLocation =\
            self.retrieveGeneLocation(geneName, self.wantedFeatureType)

        if geneLocation is None:
            print("Aborting brute force primer search: Gene name not found.")
            return

        FeatureGroup, FeaturePosition = geneLocation

        try:
            assert FeatureGroup.seq
            regionSequence = self.locateAndFetchSequence(
                FeatureGroup.seq, FeaturePosition)
        except Exception:
            print("Falling back to template search on genome file.")
            regionSequence =\
                self.locateAndFetchSequenceOnGenomeFile(
                    FeaturePosition, FeatureGroup.description)

        if not regionSequence:
            print("\n")
            print("Error: Failure on feching brute force sequence.")
            print("genomePath: %s" % self.matchedGenome)
            print("chromosome descriptor: %s" % FeatureGroup.description)
            print("location: %s" % FeaturePosition)
            return
        else:
            return regionSequence

    def launchBruteForcePrimerSearch(self, locus_name, chromosomes, Reverse):

        # BRUTE FORCE PRIMER FINDER OPERATIONS;
        geneSequenceFile = "%s.fasta" % locus_name

        PrimerSourcesDirectory = "PrimerSources"
        if not os.path.isdir(PrimerSourcesDirectory):
            os.mkdir(PrimerSourcesDirectory)

        geneSequenceFilePath =\
            os.path.join(PrimerSourcesDirectory, geneSequenceFile)

        if not os.path.isfile(geneSequenceFilePath):
            # Fetch gene sequence;
            regionSequence = self.fetchGeneSequence(locus_name)

            if regionSequence:
                # Save sequence;
                outputFile = open(geneSequenceFilePath, 'w')
                outputFile.write(str(regionSequence))
                outputFile.close()

        if os.path.isfile(geneSequenceFilePath):
            with open(geneSequenceFilePath) as f:
                geneSequenceRaw = f.read()
        else:
            print("Primer source not found.")
            return None, None

        # replace with SeqIO methods
        geneSequence = geneSequenceRaw.split("\n")
        if ">" in geneSequence[0]:
            geneSequence = geneSequence[1:]
        geneSequence = "".join(geneSequence).lower()

        foundPrimers, chr_identifier =\
            self.findPrimerBruteForce(
                chromosomes,
                geneSequence,
                Reverse=Reverse
            )

        if foundPrimers:
            print("Brute force forward primer search returns "
                  "%i primers." % len(foundPrimers))

        resultingPrimers = [p.upper() for p in foundPrimers]

        return resultingPrimers, chr_identifier

    def GeneralizePrimer(self, Primer, U):
        primer = "".join([
            "." if random.random() < U else p for p in Primer
        ])

        return primer

    def findPrimerBruteForce(self,
                             genome,
                             gene_sequence,
                             Reverse=False,
                             maximumPrimerCount=36,
                             Verbose=False):

        PRIMER_LENGTH = self.PrimerLength
        SEARCH_STEP = 5

        # FOCUS SEARCH ON A REGION ON THE MIDDLE OF THE GENE SEQUENCE;
        sequenceLength = len(gene_sequence)

        allowed_gene_sequence = gene_sequence
        if sequenceLength > self.AmpliconMaximumLength:
            HSL = sequenceLength // 2
            HSLA = self.AmpliconMaximumLength // 2
            sequenceLengthBounds = (
                HSL - HSLA,
                HSL + HSLA
            )

            allowed_gene_sequence = gene_sequence[
                sequenceLengthBounds[0]: sequenceLengthBounds[1]]

        EffectiveMinimumAmpliconLength = min(self.AmpliconMinimumLength,
                                             len(allowed_gene_sequence) // 1.2)

        FinalIndex = int(abs(
            len(allowed_gene_sequence) - EffectiveMinimumAmpliconLength
        ))

        if Reverse:
            PrimerIndexes =\
                range(len(allowed_gene_sequence) - PRIMER_LENGTH,
                      len(allowed_gene_sequence) - FinalIndex, -SEARCH_STEP)
        else:
            PrimerIndexes = range(0, FinalIndex, SEARCH_STEP)

        PrimerSequences = [
            allowed_gene_sequence[i:i + PRIMER_LENGTH]
            for i in PrimerIndexes
        ]

        if self.PrimerAllowedUncertainty:
            PrimerSequences = [
                self.GeneralizePrimer(primer, self.PrimerAllowedUncertainty)
                for primer in PrimerSequences
            ]

        # Filter primers by their real PCR capabilities if the user wants;
        if self.FindPCRViablePrimers:
            # We would evaluate the reverse complements, actually.
            PrimerReveseComplements = [
                str(Seq.Seq(p).reverse_complement())
                for p in PrimerSequences
            ]

            PrimerPCRScores = [
                (p, RealPrimers.EvaluatePrimerForPCR(rev))
                for p, rev in zip(PrimerSequences, PrimerReveseComplements)
            ]
            PrimerPCRScores = sorted(PrimerPCRScores,
                                     key=lambda ps: ps[1], reverse=True)

            PrimerSequences = [ps[0] for ps in PrimerPCRScores if ps[1] > 0]

        # Test newly aquired primers;
        foundPrimers = []
        chr_identifier = None
        for p, primer_sequence in enumerate(PrimerSequences):
            for c, _chr in enumerate(genome):
                matches, sequenceVariationName =\
                    PrimerDock.findPrimer(_chr, primer_sequence)

                if len(matches) > 1:
                    if Verbose:
                        print("Leak.")
                    continue

                if matches:
                    if chr_identifier is None:
                        chr_identifier = _chr.name
                    if Verbose:
                        print(matches[0][0].upper())
                        print(sequenceVariationName)

                    foundPrimers.append(primer_sequence)
                    if len(foundPrimers) > maximumPrimerCount:
                        return foundPrimers, chr_identifier

        return foundPrimers, chr_identifier
