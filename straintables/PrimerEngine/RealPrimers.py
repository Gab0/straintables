#!/bin/python

from Bio.SeqUtils import MeltingTemp as mt


def EvaluatePrimerForPCR(Primer):

    # -- GC content;
    gcb = [base for base in Primer.lower() if base in "gc"]
    gc = len(gcb) / len(Primer)

    if not (0.4 <= gc <= 0.6):
        return False

    # -- GC at extremities;
    extremities = Primer[:2] + Primer[-2:]
    if not all([e in "gc" for e in extremities.lower()]):
        return False

    # -- check for 2D primer formation;
    # TBD;

    # -- Melting Temperature
    MTemp = mt.Tm_NN(Primer)
    if not (52 <= MTemp <= 58):
        return False

    return True
