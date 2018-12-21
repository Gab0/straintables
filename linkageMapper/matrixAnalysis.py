# !/bin/python

import os
import pandas as pd
import numpy as np
from optparse import OptionParser

import skbio.stats.distance as skdist

parser = OptionParser()
parser.add_option("-d", dest="inputDirectory")

options, args = parser.parse_args()

# CHECK INPUT ARGUMENTS;
if not options.inputDirectory:
    print("FATAL: No input directory.")
    exit(1)

# LOAD RESULT FILES;
allFiles = os.listdir(options.inputDirectory)
arrayFiles = [File for File in allFiles if File.endswith(".aln.npy")]

arrayFilePaths = [os.path.join(options.inputDirectory, File)
                  for File in arrayFiles]

heatmaps = [np.load(filePath) for filePath in arrayFilePaths]

heatmapLabels = np.load(os.path.join(options.inputDirectory,
                                     "heatmap_labels.npy"))

Distances = [skdist.DistanceMatrix(h, heatmapLabels) for h in heatmaps]

grouping = []
for l in heatmapLabels:
    g = len(l.split("_")[-1])
    grouping.append(g)

allResults = []
for d, D in enumerate(Distances):
    print(arrayFilePaths[d])
    res = skdist.anosim(D, grouping=grouping)
    print(res)
    res["Locus"] = arrayFiles[d]
    allResults.append(res)
    print("\n")

PWM = skdist.pwmantel(Distances)
PWM_Index_Indices = [(x[0], x[1]) for x in PWM.index]
PWM_Index_Labels = [(arrayFiles[x[0]], arrayFiles[x[1]]) for x in PWM.index]
PWM.index = pd.MultiIndex.from_tuples(PWM_Index_Labels)

# INITIALIZE RESULT LISTS;
Associations = []
Mantels = []
MantelsP = []
RankingDiff = []


def matrixRankings(MATRIX):
    orderMatrix = []
    for row in MATRIX:
        row = list(row)

        order = [row.index(v) for v in sorted(row)]
        orderMatrix.append(order)
    return orderMatrix


def compareMatrixRankings(mat0, mat1):
    mat0 = matrixRankings(mat0)
    mat1 = matrixRankings(mat1)

    diff = 0
    for i in range(len(mat0)):
        for j in range(len(mat0)):
            diff += abs(mat0[i][j] - mat1[i][j])

    return diff


# ITERATE pwm_indices;
for IFA, IFB in PWM_Index_Indices:

    # FIRST EVALUATION;
    D = np.cov(heatmaps[IFA], heatmaps[IFB])
    D = sum(sum(D))
    Associations.append(D)

    # SECOND EVALUATION;
    M = skdist.mantel(heatmaps[IFA], heatmaps[IFB])

    Mantels.append(M[0])
    MantelsP.append(M[1])

    # THIRD EVALUATION;
    MDF = compareMatrixRankings(heatmaps[IFA], heatmaps[IFB])
    RankingDiff.append(MDF)

PWM["associations"] = Associations

PWM["mantel"] = Mantels
PWM["mantel_p"] = MantelsP
PWM["matrix_ranking_diff"] = RankingDiff


# cleanup dataframe;
del PWM["permutations"]

print(PWM.to_string())

for IFA, IFB in PWM_Index_Indices:
    pass

# WRITE OUTPUT PWM FILE;
pwmPath = os.path.join(options.inputDirectory, "PWMAnalysis.csv")
PWM.to_csv(pwmPath)

outputData = pd.DataFrame(allResults,
                          columns=["Locus"] + list(allResults[0].keys())[:-1])

outputPath = os.path.join(options.inputDirectory, "HeatmapAnalysis.csv")
outputData.to_csv(outputPath, index=False)

print("Wrote %s analysis file." % outputPath)
