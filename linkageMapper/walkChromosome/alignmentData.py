#!/bin/python

import os
import re
import pandas as pd
import numpy as np


class AlignmentData():
    def __init__(self, inputDirectory):
        self.dataKeys = ["Unnamed: 0", "Unnamed: 1"]
        self.loadDataFiles(inputDirectory)
        self.inputDirectory = inputDirectory

    def loadDataFile(self, inputDirectory, filename):
        FilePath = os.path.join(inputDirectory, filename)
        return pd.read_csv(FilePath)

    def loadDataFiles(self, inputDirectory):

        # LOAD RELEVANT DATABASES;
        self.PrimerData = self.loadDataFile(inputDirectory,
                                            "PrimerData.csv")

        self.PWMData = self.loadDataFile(inputDirectory,
                                         "PWMAnalysis.csv")

        self.MatchData = self.loadDataFile(inputDirectory,
                                           "MatchedRegions.csv")

        self.AlignmentData = self.loadDataFile(inputDirectory,
                                               "AlignedRegions.csv")

        # FETCH ORIGINAL HEATMAP GENOME LABELS;
        heatmapLabelsFilePath = os.path.join(
            inputDirectory,
            "heatmap_labels.npy"
        )

        self.heatmapLabels = np.load(heatmapLabelsFilePath)

        # FETCH VIEWABLE DATA INDEXES;
        OnlySequence = True
        if OnlySequence:
            last = None
            self.allowedIndexes = []
            for I in range(self.PWMData.shape[0]):
                d = self.PWMData.iloc[I]
                a = d[self.dataKeys[0]]
                if a == last:
                    continue
                self.allowedIndexes.append(I)
                last = a
        else:
            self.allowedIndexes = list(range(self.PWMData.shape[0]))

            print("Allowed: %s" % self.allowedIndexes)

    def findPWMDataRow(self, a, b):
        def setLength(w):
            return len(list(set(w)))

        for k in range(self.PWMData.shape[0]):
            d = self.PWMData.iloc[k]
            names = [d[x] for x in self.dataKeys]

            if a in names and b in names and setLength(names) == setLength([a, b]):
                print(d)
                return d

        return None

    def buildArrayPath(self, f):
        return os.path.join(self.inputDirectory, f + ".aln.npy")

    def fetchLociList(self):
        for t in self.dataKeys:
            print(self.PWMData[t])
        allLoci = [list(self.PWMData[d]) for d in self.dataKeys]
        print(allLoci)
        allLoci = [j for s in allLoci for j in s]

        return list(set(allLoci))

    def getPWMRegionIndexes(self, Index, fullName=False):
        Data = self.PWMData.iloc[Index]

        locusNames = [Data[kn] for kn in self.dataKeys]

        if not fullName:
            locusNames = [n.replace(".npy", "") for n in locusNames]

        return [self.getLocusIndex(name)
                for name in locusNames]

    def getLocusIndex(self, Name):
        Names = [
            Name,
            self.locusFromAlignmentFilename(Name)
        ]
        Results = []
        for Name in Names:
            Results += list(self.MatchData.index[self.MatchData["LocusName"] == Name])

        return Results[0]

    def locusFromAlignmentFilename(self, Filename):
        Name = re.findall("LOCI_([\w\d]+)\.", Filename)[0]
        return Name
