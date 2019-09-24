#!/bin/python

"""

straintables' main pipeline script;


"""

import os
import argparse
import pandas as pd
import shutil
import straintables
import subprocess

from Bio.Align.Applications import ClustalwCommandline

from straintables.logo import logo

from straintables.Executable import primerFinder, detectMutations,\
    compareHeatmap, matrixAnalysis

from straintables.Database import directoryManager

ClustalCommand = "clustalo"


class Options():
    def __init__(self, options):
        self.__dict__.update(options)


def find_primers(options):

    return primerFinder.Execute(options)


def run_alignment(filePrefix):
    infile = filePrefix + ".fasta"
    outfile = filePrefix + ".aln"

    aln_cmd = ClustalwCommandline(ClustalCommand, infile=infile, outfile=outfile)
    stdout, stderr = aln_cmd()

    print(stdout)

    infile = filePrefix + ".aln"
    tree_cmd = ClustalwCommandline(ClustalCommand, infile=infile, tree=True)
    tree_cmd()


def draw_tree(filePrefix):
    infile = filePrefix + ".ph"
    outfile = filePrefix + "pdf"

    treeOptions = Options({
        "InputFile": infile,
        "OutputFile": outfile
    })

    straintables.DrawGraphics.drawTree.Execute(treeOptions)


def run_meshclust(filePrefix):
    subprocess.run([
        "meshclust",
        filePrefix + ".fasta",
        "--output",
        filePrefix + ".clst",
        "--id", "0.999",
        "--align"
    ])


def detect_mutations(filePrefix):
    infile = filePrefix + ".aln"

    mutationsOptions = Options({
        "InputFile": infile,
        "PlotSubtitle": ""
    })

    detectMutations.Execute(mutationsOptions)


def matrix_analysis(WorkingDirectory):
    analysisOptions = Options({
        "WorkingDirectory": WorkingDirectory,
        "updateOnly": False
    })

    compareHeatmap.Execute(analysisOptions)
    return matrixAnalysis.Execute(analysisOptions)


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("--noamplicon", dest="DoAmplicon",
                        action="store_false", default=True)

    parser.add_argument("--noalign", dest="DoAlignment",
                        action="store_false", default=True)

    parser.add_argument("--alnmode", dest="AlignmentMode",
                        default="clustal")

    # parser.add_argument("--clustalpath", dest="ClustalPath",
    #                     default="clustalo")

    parser = primerFinder.parse_arguments(parser)
    options = parser.parse_args()

    return options


def TestMeshclust():
    # -- TEST MESHCLUST SETUP;
    if shutil.which("meshclust"):
        print("MeshClust enabled!")
        return True
    else:
        print("MeshClust not found! Disabled...")
        return False


def process_individual_region(options, locusName, MeshClustEnabled):
    filePrefix = os.path.join(options.WorkingDirectory, "LOCI_" + locusName)
    print("Running alignment for %s..." % locusName)

    run_alignment(filePrefix)
    # draw_tree(filePrefix)
    detect_mutations(filePrefix)
    if MeshClustEnabled:
        run_meshclust(filePrefix)


def Execute(options):

    if not options.PrimerFile:
        print("Fatal: No primer file specified!")
        exit(1)

    # -- SELECT WORKING DIRECTORY;
    if not options.WorkingDirectory:
        AnalysisCode = os.path.splitext(options.PrimerFile)[0]
        AnalysisCode = os.path.basename(AnalysisCode)

        WorkingDirectoryBase = "analysisResults"
        options.WorkingDirectory = os.path.join(WorkingDirectoryBase,
                                                AnalysisCode)

    # -- TEST CLUSTAL SETUP;
    if not shutil.which(ClustalCommand):
        print("%s not found! Aborting..." % ClustalCommand)
        exit(1)

    MeshClustEnabled = TestMeshclust()
    directoryManager.createDirectoryPath(options.WorkingDirectory)

    # SHOW BEAUTIFUL ASCII ART;
    print(logo)

    # -- RUN PIPELINE;
    if options.DoAmplicon:
        result = find_primers(options)
        if not result:
            print("Failure to find primers.")
            exit(1)

    AllowedAlignModes = ["clustal"]
    if options.AlignmentMode not in AllowedAlignModes:
        print("Unknown alignment mode %s." % (options.AlignmentMode))
        exit(1)

    MatchedRegions = straintables.OutputFile.MatchedRegions(options.WorkingDirectory)
    MatchedRegions.read()

    SucessfulLoci = MatchedRegions.content["LocusName"]

    if options.DoAlignment:
        for locusName in SucessfulLoci:
            process_individual_region(options, locusName, MeshClustEnabled)

    if matrix_analysis(options.WorkingDirectory):
        print("Analysis sucesfull.")


def main():
    options = parse_arguments()
    Execute(options)


if __name__ == "__main__":
    main()
