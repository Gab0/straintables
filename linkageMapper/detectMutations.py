#!/bin/python

from Bio import AlignIO
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import sys
import copy
import re
import os
import geneGraphs

import optparse


def heatmapToAxis(MATRIX, ax, labels=None):
    ax.matshow(MATRIX, cmap='gray')
    MATRIX = np.asmatrix(MATRIX)
    SIZE = MATRIX.shape[0]

    # MINOR TICKS -> GRID;
    DIV = SIZE // 3
    gridPoints = np.arange(0, SIZE, DIV)[1:-1] + 0.5

    ax.set_xticks(gridPoints, minor=True)
    ax.set_yticks(gridPoints, minor=True)

    # MAJOR TICKS -> LABELS;
    ax.set_xticks(range(SIZE))
    ax.set_yticks(range(SIZE))

    # SET LABELS;
    if labels is not None:
        ax.set_xticklabels(labels, rotation=90)
        ax.set_yticklabels(labels)


# BUILD HEATAP;
def plotHeatmap(MATRIX, sequenceNames, filename=None, subtitle=None):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    heatmapToAxis(MATRIX, ax, sequenceNames)

    # ax.grid(which='minor', color='r', linestyle='-', linewidth=2)

    if filename:
        plt.savefig(filename, bbox_inches='tight')

    watermarkLabel = re.findall("LOCI_(.+)\.aln", filename)
    # for loci similarity matrix;
    if watermarkLabel:
        watermarkLabel = watermarkLabel[0]
    # for other types of matrix;
    else:
        watermarkLabel = os.path.split(filename)[-1].split(".")[0]

    geneGraphs.watermarkAndSave(watermarkLabel, filename,
                                subtitle=subtitle, verticalLabel=340)


# reorder windows and sequenceNames;
def sortAlignments(sequenceNames, _Windows):
    def alphabeticalEngine(x):
        return x[0]

    def rflpEngine(x):
        a = x[0].split("_")[-1]
        order = len(a)
        order += len(list(set(a)))
        return order

    def stackEngine(x):
        return (rflpEngine(x), alphabeticalEngine(x))

    a = zip(sequenceNames, _Windows)
    a = sorted(a, key=stackEngine)
    sequenceNames = [x[0] for x in a]
    _Windows = [x[1] for x in a]
    return sequenceNames, _Windows


def storeMatrixData(MATRIX, baseFileName, subtitle=None):
    # SAVE HEATMAP GRAPHIC;
    plotHeatmap(MATRIX, sequenceNames, baseFileName + ".pdf", subtitle=subtitle)

    # SAVE HEATMAP MATRIX;
    np.save(baseFileName, MATRIX)

    labelsPath = os.path.join(os.path.dirname(baseFileName), "heatmap_labels")
    np.save(labelsPath, np.array(sequenceNames))


if __name__ == "__main__":

    parser = optparse.OptionParser()

    parser.add_option('-i', dest="InputFile")
    parser.add_option('-s', dest="PlotSubtitle")

    options, args = parser.parse_args()

    alignPath = options.InputFile
    if not alignPath:
        print("No input file specified!")
        exit(1)

    # LOAD ALIGNMENT;
    Alignment = AlignIO.read(open(alignPath), "clustal")

    # GET SEQUENCE NAMES;
    sequenceNames = [d.id for d in Alignment]

    # BUILD VARIATION WINDOWS;
    Windows = []
    InfoWindows = []
    for s in range(len(Alignment[0].seq)):
        window = [dna.seq[s] for dna in Alignment]
        if len(list(set(window))) > 1:
            print(window)
            Windows.append(window)
            InfoWindows.append([s] + window)

            snp_variations = len(set(window))
            if snp_variations > 2:
                print("TRI_SNP")

    # VARIATION WINDOWS TO SIMILARITY MATRIX;
    # transpose windows;
    _Windows = list(zip(*Windows))

    sequenceNames, _Windows = sortAlignments(sequenceNames, _Windows)

    mat = range(len(_Windows))
    nb_snp = len(_Windows[0])
    print(len(Windows))
    print(nb_snp)
    MATRIX = [[0 for j in mat] for i in mat]
    for i in mat:
        for j in mat:
            similarity = 0
            for k in range(nb_snp):
                A = _Windows[i][k]
                B = _Windows[j][k]
                if A == B:
                    similarity += 1
            similarity /= nb_snp
            MATRIX[i][j] = similarity

    d = np.matrix(MATRIX)

    # SHOW SIMILARITY MATRIX;
    print(d.shape)
    for i in MATRIX:
        for j in i:
            print("%.2f" % j, end=' ')
        print("\n")

    # PROCESS MATRIX;
    c = np.exp(d)
    _MATRIX = copy.deepcopy(MATRIX)
    for i in mat:
        for j in mat:
            v = MATRIX[i][j]
            new_v = (v / np.sum(MATRIX[i]))
            _MATRIX[i][j] = new_v

    # SAVE DIVERSE VERSIONS OF THE MATRIX;
    storeMatrixData(MATRIX, alignPath, subtitle=options.PlotSubtitle)
    #storeMatrixData(_MATRIX, alignPath + "alt")

    # SAVE SNP INFO DATA FILE;
    DATA = pd.DataFrame(InfoWindows, columns=["POS"] + sequenceNames)
    DATA.to_csv(alignPath + ".csv", index=False)
