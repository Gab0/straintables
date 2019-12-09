#!/bin/python

import os
from Bio import SeqIO


class Gene():
    def __init__(self, Name, hasName):
        self.Name = Name
        self.hasName = hasName


def loadGenesFromScaffoldAnnotation(Scaffold):
    genes = []
    for feature in Scaffold.features:
        COND1 = "gene" in feature.qualifiers.keys()
        COND2 = "locus_tag" in feature.qualifiers.keys()

        NAME = None
        hasName = False

        if COND1:
            NAME = feature.qualifiers["gene"][0]
        elif COND2:
            NAME = feature.qualifiers["locus_tag"][0]

        # NAME or COND1?
        if NAME:
            genes.append(Gene(NAME, hasName))

    return genes


def loadFeatures(annotationFilePath):
    annotation = SeqIO.read(annotationFilePath, "genbank")
    outputFeatures = []
    for feature in annotation.features:
        if feature.type in ["gene"]:
            outputFeatures.append(feature)

    return outputFeatures


def filterScaffoldsByIdentifier(annotationScaffolds, identifier):
    wantedIdentifiers = [
        identifier,
        "chromosome_%s" % identifier
    ]

    # -- pick only those that contain the identifier
    allIdentifiers = []
    wantedScaffolds = []
    for Scaffold in annotationScaffolds:
        Qualifiers = Scaffold.features[0].qualifiers
        if 'chromosome' in Qualifiers.keys():
            ChromosomeName = Qualifiers['chromosome'][0]
            allIdentifiers.append(ChromosomeName)
            for wantedIdentifier in wantedIdentifiers:
                if wantedIdentifier.lower() == ChromosomeName.lower():
                    wantedScaffolds.append(Scaffold)

    return wantedScaffolds


def loadAnnotation(annotationFolder, identifier=None, Verbose=False):

    """

    Returns the most suitable annotation file from the annotations folder,
    as a list of scaffolds.

    """

    annotationFiles = os.listdir(annotationFolder)
    annotationFiles = sorted([File
                              for File in annotationFiles
                              if File.endswith(".gbff")])
    annotationFilePaths = [
        os.path.join(annotationFolder, annotationFile)
        for annotationFile in annotationFiles
    ]

    if not annotationFiles:
        print("Annotation file not found! Check your annotation folder.")
        exit(1)

    def sortScaffold(scaffold):
        genes = loadGenesFromScaffoldAnnotation(scaffold)
        return len(genes)

    def scoreAnnotation(annotation):
        genes = loadGenesFromScaffoldAnnotation(annotation[0])
        return len(genes)

    best_score = 0
    chosen = (None, None)
    for annotationFilePath in annotationFilePaths:
        annotationScaffolds = list(SeqIO.parse(annotationFilePath, "genbank"))

        if identifier:
            annotationScaffolds =\
                filterScaffoldsByIdentifier(
                    annotationScaffolds,
                    identifier
                )

        if annotationScaffolds:
            annotationScaffolds = sorted(
                annotationScaffolds,
                key=sortScaffold,
                reverse=True
            )

            score = scoreAnnotation(annotationScaffolds)
            if score > best_score:
                chosen = (annotationFilePath, annotationScaffolds)
                best_score = score
                if Verbose:
                    print(annotationFilePath)
                    print(score)

    if Verbose:
        print("\n====")
        for aS in annotationScaffolds:
            print(len(aS.features))

    if chosen == (None, None):
        print("Warning: No suitable annotation file found.")

    return chosen
