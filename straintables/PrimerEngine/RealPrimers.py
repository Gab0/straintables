#!/bin/python

"""

This module deals with evaluating the
usability of primer sequences on real world PCR,
calculating parameters such as Tm and GC content.


"""
from typing import List
from collections import OrderedDict
from Bio.SeqUtils import MeltingTemp as mt
import primer3


def ValueInBounds(Value: float, Target: float,
                  Tolerance: float, exp: float = 0.4):
    return max(abs(Target - Value) - Tolerance, 0) ** exp


def PenaltyGCContent(Primer):

    # -- GC content;
    gcb = [base for base in Primer.lower() if base in "gc"]
    gc = len(gcb) / len(Primer)

    TargetGC = 0.5
    ToleranceGC = 0.1

    Penalty = ValueInBounds(gc, TargetGC, ToleranceGC)

    return Penalty


def PenaltyGCExtremities(Primer):

    # -- GC at extremities;
    extremities = Primer[:2] + Primer[-2:]
    scr = sum([e in "gc" for e in extremities.lower()]) / len(extremities)
    Penalty = 1.0 - scr
    return Penalty


def PenaltyMeltingTemperature(Primer: str):

    # -- Melting Temperature
    MTemp = mt.Tm_NN(Primer)

    Penalty = ValueInBounds(MTemp, 55, 3)

    return Penalty


def BuildPrimerReport(PrimerIdentifier: str, PrimerPair: List[str]):
    IDs = ["5'", "3'"]

    def createPrimer(INP):
        (i, Primer) = INP
        return OrderedDict({
            "Primer": "%s %s" % (PrimerIdentifier, IDs[i]),
            "Sequence": PrimerPair[i],
            "GCContent": PenaltyGCContent(Primer),
            "GCExtremities": PenaltyGCExtremities(Primer),
            "Tm p": PenaltyMeltingTemperature(Primer),
            "Tm": mt.Tm_NN(Primer)
        })

    PrimerReport = list(map(createPrimer, enumerate(PrimerPair)))
    deltaT = abs(PrimerReport[0]["Tm"] - PrimerReport[1]["Tm"])
    deltaTpenalty = ValueInBounds(deltaT, 0, 5)

    for k, _ in enumerate(PrimerReport):
        PrimerReport[k]["Tm p"] -= deltaTpenalty

    return PrimerReport


def CalculatePrimerPairScore(PrimerReport: List[OrderedDict]):
    Score = 1.0

    penaltyKeys = ["GCContent", "GCExtremities", "Tm p"]

    for pK in penaltyKeys:
        Score -= sum([d[pK] for d in PrimerReport])

    return Score


def EvaluatePrimerForPCR(Primer: str):
    """

    Legacy method of calculating primer PCR score;

    """

    Score = 1.0

    Score -= PenaltyGCContent(Primer)
    Score -= PenaltyGCExtremities(Primer)
    Score -= PenaltyMeltingTemperature(Primer)

    # -- check for 2D primer formation;
    Hairpin = primer3.calcHairpin(Primer)
    if Hairpin.structure_found:
        Score -= 1.0

    return Score
