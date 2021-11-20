#!/bin/python

import numpy as np


def buildVariationWindows(Alignment, Verbose=False):
    Windows = []
    InfoWindows = []

    for s in range(len(Alignment[0].seq)):
        letter_set = [dna.seq[s] for dna in Alignment]
        if len(list(set(letter_set))) > 1:

            # Special Cases: HUGE INSERTIONS
            # AT THE BEGINNING OR END OF ANY SEQUENCE;
            # TBD
            # print(letter_set)
            Windows.append(letter_set)
            InfoWindows.append([s] + letter_set)

            snp_variations = len(set(letter_set))
            if snp_variations > 2:
                if Verbose:
                    print("TRIALLELIC_SNP")

    return Windows, InfoWindows


def normalize_distances(Matrix):
    """Smart normalization to improve Matrix visualization.

    :param Matrix: Input matrix.
    :returns: Modified matrix.

    """
    original_value = Matrix[0][0]
    for k in np.nditer(Matrix):
        if k not in [0, 1]:
            W = k
            break

    np.fill_diagonal(Matrix, W)
    xmax, xmin = Matrix.max(), Matrix.min()
    Matrix = (Matrix - xmin)/(xmax - xmin)

    np.fill_diagonal(Matrix, original_value)


def buildMatrixFromWindow(Alignment, _Windows, Verbose=False):
    def isLetter(v, includeN=True):
        v = v.lower()
        L = ['a', 'c', 't', 'g', 'u']

        if includeN:
            L += ['n']

        if v in L:
            return True
        return False

    def isUnknown(v):
        v = v.lower()
        if v == 'n':
            return True
        return False

    def isGap(v):
        if v == '-':
            return True
        return False

    # PROCESS VARIATION WINDOWS;
    n = len(_Windows)
    a = len(Alignment)
    mat = range(n)

    if Verbose:
        print(n)
        print(_Windows)

    MATRIX = np.matrix(np.zeros((n, n)))

    if not n:
        return np.matrix(np.zeros((a, a))), 0

    nb_possible_snp = len(_Windows[0])

    if Verbose:
        print(nb_possible_snp)

    # CALCULATE TOTAL NUMBER OF SNPS;
    snp_idx = []
    for k in range(nb_possible_snp):
        c = set([window[k] for window in _Windows])
        if sum(map(lambda x: isLetter(x, includeN=False), c)) > 1:
            snp_idx.append(k)
    nb_snp = len(snp_idx)

    for i in mat:
        for j in mat:
            # if i > j:
            #    MATRIX[i, j] = MATRIX[j, i]
            #    continue
            if j == i:
                MATRIX[i, j] = 1
                continue

            # CALCULATE SIMILARITY;
            similarity = 0
            for k in snp_idx:
                A = _Windows[i][k]
                B = _Windows[j][k]
                if any(
                        not isLetter(x) and not isGap(x)
                        for x in (A, B)):
                    continue

                if A.lower() == B.lower() or \
                   any(isUnknown(x) for x in (A, B)):
                    similarity += 1
            try:
                # FIXME: Which similarity function to use?
                similarity = similarity / len(Alignment[0])
            except ZeroDivisionError:
                print("ZDO")
                similarity = 1
            MATRIX[i, j] = similarity
    # Convert similarity to dissimilarity;
    MATRIX = 1 - MATRIX

    return MATRIX, nb_snp
