#!/bin/python
import os
from argparse import ArgumentParser

import straintables
from straintables.Executable import GenomePipeline

from Bio import SeqIO

Description = """

Secondary straintables' pipeline.
Will initialize an output directory based on
a multifasta sequence file per wanted region.

The resulting files can be inspected manually
or used to build dissimilarity matrices through
the 'stview' command.

"""


def Execute(options):
    allFiles = os.listdir(options.WorkingDirectory)

    MeshClustEnabled = GenomePipeline.TestMeshclust()

    AllRegions = straintables.OutputFile.MatchedRegions(
        options.WorkingDirectory)

    AllRegionsData = []
    FakePrimerData = []
    for File in allFiles:
        if not File.endswith(".fasta"):
            continue
        if File.startswith(straintables.Definitions.FastaRegionPrefix):
            continue

        filePrefix = os.path.splitext(File)[0]
        filePath = os.path.join(options.WorkingDirectory, File)

        outputFilePath = os.path.join(
            options.WorkingDirectory,
            "%s%s" % (straintables.Definitions.FastaRegionPrefix, File))

        sequences = SeqIO.parse(filePath, format="fasta")

        SeqIO.write(sequences, open(outputFilePath, 'w'), format="fasta")

        GenomePipeline.process_individual_region(options,
                                           filePrefix,
                                           MeshClustEnabled)

        FakePrimerData.append({
            "Locus": filePrefix,
            "Sequence": "NOSEQ-FROMFASTA",
            "PositionStart": 0,
            "PositionEnd": 0
        })
        AllRegionsData.append({
            "LocusName": filePrefix,
            "RebootCount": 0,
            "MeanLength": 0,
            "StartPosition": 0
        })

    PrimerData = straintables.OutputFile.PrimerData(options.WorkingDirectory)
    PrimerData.add(FakePrimerData)
    PrimerData.write()

    AllRegions.add(AllRegionsData)
    AllRegions.write()

    if GenomePipeline.matrix_analysis(options.WorkingDirectory):
        print("Sucess.")


def main():
    parser = ArgumentParser(description=Description)
    parser.add_argument("-d", dest="WorkingDirectory")

    parser.add_argument("--clustalpath", dest="ClustalPath",
                        default=straintables.Definitions.ClustalCommand)
    options = parser.parse_args()

    Execute(options)


if __name__ == "__main__":
    main()
