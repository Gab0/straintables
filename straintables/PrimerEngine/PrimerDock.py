#!/bin/python

import random
import re
import itertools
import numpy as np
import copy
from Bio.Seq import Seq

from . import GeneticEntities


PrimerTypes = ["ForwardPrimer", "ReversePrimer"]


"""

findPrimer(string, string):
returns the list of matched primers along with the sequence modification
 on the primer that resulted in match.

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


def searchPrimerPairOnGenome(locusName, primerPair, genome):

    # INIT SUPPORT LIST;
    matchedPrimers = {}

    matchSuccess = [True, True]

    # ITERATE THRU PRIMER TYPES;
    for PT, PrimerType in enumerate(PrimerTypes):

        # screen info purposes;
        print(">> Matching %s" % PrimerType)

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
            # outright delete matches that
            # doesn't have a pair in the same chromosome;
            # if PrimerTypes[1 - PT] in matchedPrimers.keys():
            opponentIndexes = [
                P.chr_index
                for P in matchedPrimers[PrimerTypes[1 - PT]]
            ]
            for p, P in enumerate(matchedPrimers[PrimerType]):
                if P.chr_index not in opponentIndexes:
                    matchedPrimers[PrimerType][p] = None

            matchedPrimers[PrimerType] = [
                p for p in matchedPrimers[PrimerType] if p
            ]

            # might be unnecessary
            matchedPrimers[PrimerType] = sorted(matchedPrimers[PrimerType],
                                                key=lambda x: x.chr_length,
                                                reverse=True)

            # PRIMER LEAK IS UNAVOIDABLE;
            # !! This len is always 1 at this point (read above);
            # if len(matchedPrimers[PrimerType]) > 1:
            #    matchSuccess[PT] = False
        print()

    print("\n\n")

    # RETRIEVE INTERPRIMER SEQUENCES;

    # MATCHES ON THE SAME CHROMOSOME.. IDEAL SCENARIO;
    mp = (matchedPrimers[PrimerTypes[0]], matchedPrimers[PrimerTypes[1]])
    if all(mp):
        mp = [p[0] for p in mp]
        if mp[0].chr_length == mp[1].chr_length:

            # if anything goes wrong while building the amplicon,
            # the match fails.
            # Amplicon may rase errors deliberately
            # when match result is not optimal.
            try:
                amplicon = GeneticEntities.Amplicon(genome, *mp)
            except ValueError:
                return "", [False, False]

            return amplicon.Sequence, mp

    return "", matchSuccess


def matchPrimerOnGenome(genome, PrimerSequence,
                        PrimerType, maxNumberMatches=None):
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

            # SEARCH REGION ON ANNOTATED CHROMOSOMES;
            print("Found at chromosome %s" % chromosome.name)
            print("pos %i of %i" % (match.start(), chromosome.length))

            # RECORD PRIMER;
            matchedPrimer = GeneticEntities.primerMatch(
                match,
                PrimerType,
                chromosome,
                PrimerSequence
            )

            matchedPrimers.append(matchedPrimer)

            if maxNumberMatches is not None:
                if len(matchedPrimer) >= maxNumberMatches:
                    break

    return matchedPrimers


def validatePrimer(V):
    if type(V) == np.float64:
        return False
    elif type(V) == float:
        return False
    elif not V:
        return False
    else:
        return True


def matchLocusOnGenomes(locus_name,
                        locus_info,
                        genomes,
                        overallProgress=(0, 0),
                        rebootTolerance=20,
                        allowN=False,
                        bruteForceSearcher=None):
    LocusAmpliconSet = {}

    def initPrimerQueuer():
        return {k: [] for k in locus_info.keys()}

    primerTrash = initPrimerQueuer()
    testablePrimers = initPrimerQueuer()

    # load primer pair data from user-made Primer file;
    primerPair = dict(locus_info)

    RebootCount = 0
    chr_identifier = None

    FailingGenomes = [Genome.name for Genome in genomes]

    assert(len(list(set(FailingGenomes))) == len(FailingGenomes))
    # ITERATE GENOMES UNTIL SUCCESS;
    for Genome in itertools.cycle(genomes):
        # print("Genome: %s --\n" % genome)

        # Show header;
        HEADER_INFO = (
            *overallProgress,
            RebootCount + 1,
            locus_name,
            Genome.name
        )

        Message = ">>> Region %i of %i | run number %i ->  Searching %s on %s"
        print(Message % HEADER_INFO)

        # primer sequence may be uknown/invalid;
        primerIntegrity = [
            validatePrimer(primerPair[PrimerType])
            for PrimerType in PrimerTypes
        ]

        AmpliconSequence = ""
        print("\n")

        # fix this later..
        MatchedPrimers = copy.deepcopy(primerIntegrity)

        if all(primerIntegrity):
            print("Searching sequence for locus %s" % locus_name)
            AmpliconSequence, MatchedPrimers =\
                searchPrimerPairOnGenome(locus_name, primerPair, Genome)

            # GET CHR NAME ONLY FOR THE FIRST GENOME;
            if Genome.name == genomes[0].name:
                if type(MatchedPrimers[0]) == GeneticEntities.primerMatch:
                    chr_identifier = MatchedPrimers[0].chr_name

        # FAILED TO MATCH A PRIMER?
        if not all(MatchedPrimers):
            for PT, match in enumerate(MatchedPrimers):
                PrimerType = PrimerTypes[PT]
                if not match:
                    print("Resetting %s!" % PrimerType)

                    # delete current primer;
                    if validatePrimer(primerPair[PrimerType]):
                        primerTrash[PrimerType].append(primerPair[PrimerType])

                    primerPair[PrimerType] = None

                    # RESET ALL MATCHES!
                    LocusAmpliconSet = {}

                    RebootCount += 1
                    # brute force for the problematic genome!
                    print("Searching new primer on gene sequence...")
                    # print(bruteForceSearcher.matchedGenome)

                    # MANAGE PRIMERS ON QUEUE;
                    if not testablePrimers[PrimerType]:
                        if bruteForceSearcher:
                            BF = bruteForceSearcher.launchBruteForcePrimerSearch(locus_name, Genome, PT)
                            newPrimers, _chr_identifier = BF
                            if not newPrimers:
                                print("Warning: No bruteforce primers found...")

                            testablePrimers[PrimerType] = newPrimers
                        else:
                            FailureReason = "Brute force primer finder is disabled."
                            print(FailureReason)
                            return RegionMatchFailure(FailureReason)

                    while testablePrimers[PrimerType]:
                        print("Fetching new primer from reserve...")
                        random.shuffle(testablePrimers[PrimerType])
                        primerPair[PrimerType] = testablePrimers[PrimerType].pop()
                        if primerPair[PrimerType] not in primerTrash[PrimerType]:
                            break
        else:
            # Final sucess check.
            Sucess = True
            print("Found amplicon of length %i" % len(AmpliconSequence))

            if len(AmpliconSequence) <= 100:
                # AMPLICON IS TOO SHORT;
                print("Resetting due to short amplicon.")
                Sucess = False

            if allowN and "n" in AmpliconSequence.lower():
                print("Resetting due to N found in sequence.")
                Sucess = False

            if Sucess:
                LocusAmpliconSet[Genome.name] = AmpliconSequence
                if Genome.name in FailingGenomes:
                    FailingGenomes.remove(Genome.name)
            else:
                print("Resetting all primers, amplicon too short.")
                for s in range(2):
                    PP = primerPair[PrimerTypes[s]]
                    if validatePrimer(PP):
                        primerTrash[PrimerTypes[s]].append(PP)
                    primerPair[PrimerTypes[s]] = None
                # Reset Matches;
                LocusAmpliconSet = {}

        # CHECK LOCUS PROGRESS ACROSS ALL GENOMES;
        progressMark = len(list(LocusAmpliconSet.keys())) / len(genomes)
        print("%s match progress: %.2f%%" % (locus_name, progressMark * 100))
        print()

        # IS REGION DONE?
        if progressMark == 1:
            print(">>> Matched %s on all genomes." % locus_name)
            print()
            print()

            primerPair["RebootCount"] = RebootCount

            # exit loop...
            return RegionMatchSuccess(LocusAmpliconSet,
                                      MatchedPrimers, chr_identifier)

        # IF IT FAILS ENOUGH, SKIP LOCUS;
        elif RebootCount > rebootTolerance:
            FailureReason = "Too many reboots."
            return RegionMatchFailure(FailureReason,
                                      FailedGenomes=FailingGenomes)


class RegionMatchFailure():
    def __init__(self, reason, FailedGenomes=None):
        self.reason = reason
        self.FailedGenomes = FailedGenomes


class RegionMatchSuccess():
    def __init__(self,
                 LocusAmpliconSet=None,
                 MatchedPrimers=None,
                 chr_identifier=None):
        self.LocusAmpliconSet = LocusAmpliconSet
        self.MatchedPrimers = MatchedPrimers
        self.chr_identifier = chr_identifier
