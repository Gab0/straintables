#!/bin/python

import numpy as np
import sys
import os

from optparse import OptionParser

sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))

from . import detectMutations


def Execute(options):
    allFiles = os.listdir(options.WorkingDirectory)
    arrayFilePaths = [os.path.join(options.WorkingDirectory, File)
                      for File in allFiles if File.endswith(".aln.npy")]

    heatmaps = [np.load(filePath) for filePath in arrayFilePaths]

    heatmapLabels = np.load(os.path.join(options.WorkingDirectory,
                                         "heatmap_labels.npy"))

    heatmap = 1 - np.abs(heatmaps[0] - heatmaps[1])

    outputPath = os.path.join(options.WorkingDirectory, "discrepancy_matrix.pdf")
    detectMutations.createPdfHeatmap(heatmap, heatmapLabels, outputPath)


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-d", dest="WorkingDirectory")

    options, args = parser.parse_args()

    Execute(options)
