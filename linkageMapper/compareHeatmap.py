#!/bin/python

import numpy as np
import sys
import os

import detectMutations

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", dest="inputDirectory")

options, args = parser.parse_args()


allFiles = os.listdir(options.inputDirectory)
arrayFilePaths = [os.path.join(options.inputDirectory, File)
                  for File in allFiles if File.endswith(".aln.npy")]

heatmaps = [np.load(filePath) for filePath in arrayFilePaths]


heatmapLabels = np.load(os.path.join(options.inputDirectory, "heatmap_labels.npy"))
heatmap = 1 - np.abs(heatmaps[0] - heatmaps[1])

outputDir = os.path.dirname(sys.argv[1])
outputPath = os.path.join(outputDir, "discrepancy_matrix.pdf")
detectMutations.plotHeatmap(heatmap, heatmapLabels, outputPath)
