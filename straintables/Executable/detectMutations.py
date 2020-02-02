#!/bin/python

from Bio import AlignIO
import pandas as pd
import numpy as np

import shutil
import copy
import os

from straintables import DistanceEngine, Definitions
import optparse


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


def storeMatrixData(MATRIX, baseFileName, sequenceNames, subtitle=None):
    # SAVE HEATMAP MATRIX;
    np.save(baseFileName, MATRIX)

    labelsPath = os.path.join(os.path.dirname(baseFileName), "heatmap_labels")
    np.save(labelsPath, np.array(sequenceNames))


def updateAlignmentInfo(AlignmentInfoFilepath, data):
    AlignmentInfoColumns = [
            "LocusName",
            "AlignmentLength",
            "SNPCount"
        ]

    AlignmentInfo = None
    if os.path.isfile(AlignmentInfoFilepath):
        AlignmentInfo = pd.read_csv(AlignmentInfoFilepath)
        if list(AlignmentInfo.columns) != AlignmentInfoColumns:
            AlignmentInfo = None

    if AlignmentInfo is None:
        AlignmentInfo = pd.DataFrame(columns=AlignmentInfoColumns)

    if "LocusName" in AlignmentInfo.keys():
        ExistingAlignmentEntry = AlignmentInfo[
            AlignmentInfo["LocusName"] == data["LocusName"]
        ]
        if not ExistingAlignmentEntry.empty:
            pass

    AlignmentInfo = AlignmentInfo.append([data], ignore_index=True)
    AlignmentInfo.to_csv(AlignmentInfoFilepath, index=False)


def Execute(options):

    alignPath = options.InputFile

    if not alignPath:
        print("No input file specified!")
        exit(1)

    # LOAD ALIGNMENT;
    Alignment = AlignIO.read(open(alignPath), "clustal")

    # GET SEQUENCE NAMES;
    sequenceNames = [d.id for d in Alignment]

    # BUILD VARIATION WINDOWS;
    Windows, InfoWindows = DistanceEngine.buildVariationWindows(Alignment)

    # VARIATION WINDOWS TO SIMILARITY MATRIX;
    # transpose windows;
    _Windows = list(zip(*Windows))

    # if _Windows is empty (no variation on group),
    # do not reorder those sequenceNames.
    if _Windows:
        sequenceNames, _Windows = sortAlignments(sequenceNames, _Windows)

    MATRIX, nb_snp = DistanceEngine.buildMatrixFromWindow(Alignment, _Windows)

    # CHECK MATRIX HEALTH
    globalMean = np.mean(MATRIX)
    healthFail = False
    matrixRowMeans = []
    for row in MATRIX:
        rm = np.mean(row)
        matrixRowMeans.append(rm)
        if rm < (globalMean / 2):
            healthFail = True

    mat = range(len(_Windows))
    # SHOW SIMILARITY MATRIX;
    print(MATRIX.shape)
    for i in mat:
        for j in mat:
            print("%.2f" % MATRIX[i, j], end=' ')
        print("\n")

    # PROCESS MATRIX;
    _MATRIX = copy.deepcopy(MATRIX)
    for i in mat:
        for j in mat:
            v = MATRIX[i, j]
            new_v = (v / np.sum(MATRIX[i]))
            _MATRIX[i, j] = new_v

    # -- Write Matrix file;
    storeMatrixData(MATRIX,
                    alignPath,
                    sequenceNames,
                    subtitle=options.PlotSubtitle)

    # SAVE SNP INFO DATA FILE;
    DATA = pd.DataFrame(InfoWindows, columns=["POS"] + sequenceNames)
    SnpInfoPath = alignPath + ".csv"

    DATA.to_csv(SnpInfoPath, index=False)
    print("Written mutation region info at %s." % SnpInfoPath)

    FileDirectory = os.path.dirname(options.InputFile)
    AlignmentInfoFilepath = os.path.join(FileDirectory, "AlignedRegions.csv")
    LocusName = os.path.splitext(os.path.split(alignPath)[-1])[0].replace(
        Definitions.FastaRegionPrefix, "")

    # Update Alignment Info database;
    data = {
        "LocusName": LocusName,
        "AlignmentLength": len(Alignment[0]),
        "SNPCount": nb_snp
    }
    updateAlignmentInfo(AlignmentInfoFilepath, data)

    # SHOW OUTCOMES
    if healthFail:
        print("SIMILARITY MATRIX IS NOT HEALTHY.")
        problematicDirectory = os.path.join(FileDirectory, "problematic")
        if not os.path.isdir(problematicDirectory):
            os.mkdir(problematicDirectory)
        shutil.copyfile(options.InputFile,
                        os.path.join(problematicDirectory,
                                     os.path.basename(options.InputFile)))
    print(matrixRowMeans)


if __name__ == "__main__":

    parser = optparse.OptionParser()

    parser.add_option('-i', dest="InputFile")
    parser.add_option('-s', dest="PlotSubtitle")

    options, args = parser.parse_args()

    Execute(options)
