#!/bin/python
import os
from argparse import ArgumentParser

import straintables
from straintables.Executable import Pipeline


def Execute(options):
    allFiles = os.listdir(options.WorkingDirectory)

    MeshClustEnabled = Pipeline.TestMeshclust()

    AllRegions = straintables.OutputFile.MatchedRegions(
        options.WorkingDirectory)

    AllRegionsData = []
    for File in allFiles:
        if not File.endswith(".fasta"):
            continue
        filePrefix = os.path.splitext(File)[0]

        Pipeline.process_individual_region(options,
                                           filePrefix,
                                           MeshClustEnabled)

        AllRegionsData.append({
            "LocusName": filePrefix,
            "RebootCount": 0,
            "MeanLength": 0,
            "StartPosition": 0
        })

    AllRegions.add(AllRegionsData)
    AllRegions.write()

    if Pipeline.matrix_analysis(options.WorkingDirectory):
        print("Sucess.")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("-d", dest="WorkingDirectory")

    options, args = parser.parse_args()

    Execute(options)
