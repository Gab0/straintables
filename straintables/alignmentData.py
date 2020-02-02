#!/bin/python

import os
import re
import numpy as np

from . import OutputFile, Definitions


class AlignmentData():
    def __init__(self, inputDirectory):
        self.dataKeys = ["Unnamed: 0", "Unnamed: 1"]
        self.loadDataFiles(inputDirectory)
        self.inputDirectory = inputDirectory

    def loadDataFile(self, inputDirectory, FileType):
        data = FileType(inputDirectory)
        data.read()
        return data.content

    def loadDataFiles(self, inputDirectory):

        # LOAD RELEVANT DATABASES;
        self.PrimerData = self.loadDataFile(inputDirectory,
                                            OutputFile.PrimerData)

        self.PWMData = self.loadDataFile(inputDirectory,
                                         OutputFile.PWMAnalysis)

        self.MatchData = self.loadDataFile(inputDirectory,
                                           OutputFile.MatchedRegions)

        self.AlignmentData = self.loadDataFile(inputDirectory,
                                               OutputFile.AlignedRegions)

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
        # else:
        #    self.allowedIndexes = list(range(self.PWMData.shape[0]))
        #    print("Allowed: %s" % self.allowedIndexes)

    def findPWMDataRow(self, a_name, b_name):
        def setLength(w):
            return len(list(set(w)))

        for k in range(self.PWMData.shape[0]):
            d = self.PWMData.iloc[k]
            names = [d[x] for x in self.dataKeys]

            fullname = "".join(names)
            if a_name != b_name:
                if a_name in fullname:
                    if b_name in fullname:
                        if setLength(names) == setLength([a_name, b_name]):
                            print(d)
                            return d

        return None

    def buildArrayPath(self, f):
        possibleFilenames = [
            "%s.aln.npy" % f,
            "%s%s.aln.npy" % (Definitions.FastaRegionPrefix, f)
        ]

        possibleFilepaths = [
            os.path.join(self.inputDirectory, f)
            for f in possibleFilenames
        ]

        for filepath in possibleFilepaths:
            if os.path.isfile(filepath):
                return filepath

        print("Failure to find array %s" % f)
        exit(1)

    def fetchOriginalLociList(self):
        return list(self.AlignmentData["LocusName"])

    def fetchLociList(self):
        # for t in self.dataKeys:
        #    print(self.PWMData[t])
        allLoci = [list(self.PWMData[d]) for d in self.dataKeys]
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
            D = self.MatchData.index[self.MatchData["LocusName"] == Name]
            Results += list(D)

        return Results[0]

    def getRegionNamesFromIndex(self, regionIndexes):
        return [
            self.MatchData["LocusName"].iloc[i]
            for i in regionIndexes
        ]

    def getMatchDataFromNames(self, regionNames):
        MatchData = [
            self.MatchData[
                self.MatchData.LocusName == name
            ]
            for name in regionNames
        ]

        return [md.iloc[0] for md in MatchData if not md.empty]

    def locusFromAlignmentFilename(self, Filename):
        Name = re.findall(r"%s([\w\d]+)\." % Definitions.FastaRegionPrefix,
                          Filename)[0]
        return Name
