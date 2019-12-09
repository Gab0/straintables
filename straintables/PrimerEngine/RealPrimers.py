#!/bin/python

from Bio.SeqUtils import MeltingTemp as mt


def ValueInBounds(Value, Target, Tolerance, exp=0.4):
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


def PenaltyMeltingTemperature(Primer):

    # -- Melting Temperature
    MTemp = mt.Tm_NN(Primer)

    Penalty = ValueInBounds(MTemp, 55, 3)

    return Penalty


def EvaluatePrimerForPCR(Primer):

    Score = 1.0

    Score -= PenaltyGCContent(Primer)
    Score -= PenaltyGCExtremities(Primer)
    Score -= PenaltyMeltingTemperature(Primer)

    # -- check for 2D primer formation;
    # TBD;

    return Score
