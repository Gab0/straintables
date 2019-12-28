#!/bin/python

import numpy as np


def buildVariationWindows(Alignment):
    Windows = []
    InfoWindows = []

    for s in range(len(Alignment[0].seq)):
        letter_set = [dna.seq[s] for dna in Alignment]
        if len(list(set(letter_set))) > 1:

            # Special Cases: HUGE INSERTIONS
            # AT THE BEGINNING OR END OF ANY SEQUENCE;
            # TBD
            print(letter_set)
            Windows.append(letter_set)
            InfoWindows.append([s] + letter_set)

            snp_variations = len(set(letter_set))
            if snp_variations > 2:
                print("TRIALLELIC_SNP")

    return Windows, InfoWindows


def buildMatrixFromWindow(Alignment, _Windows):
    def isLetter(v):
        v = v.lower()
        if v in ['a', 'c', 't', 'g', 'u', 'n']:
            return True
        return False

    def isUknown(v):
        v = v.lower()
        if v == 'n':
            return True
        return False

    # PROCESS VARIATION WINDOWS;
    mat = range(len(_Windows))

    print(len(_Windows))
    print(_Windows)

    MATRIX = [[0 for j in mat] for i in mat]

    if MATRIX:
        nb_snp = len(_Windows[0])
        print(nb_snp)
        for i in mat:
            for j in mat:
                similarity = 0
                for k in range(nb_snp):
                    A = _Windows[i][k]
                    B = _Windows[j][k]

                    if A == B:
                        similarity += 1
                    elif isUknown(A) and isLetter(B):
                        similarity += 1
                    elif isUknown(B) and isLetter(A):
                        similarity += 1

                similarity = similarity / nb_snp

                MATRIX[i][j] = 1 - similarity

        MATRIX = np.matrix(MATRIX)

    else:
        MATRIX = np.matrix(np.zeros((len(Alignment), len(Alignment))))
        nb_snp = 0

    return MATRIX, nb_snp
