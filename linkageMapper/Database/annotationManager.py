#!/bin/python

import os
from Bio import SeqIO


def loadFeatures(annotationFilePath):
    annotation = SeqIO.read(annotationFilePath, "genbank")
    outputFeatures = []
    for feature in annotation.features:
        if feature.type in ["gene"]:
            outputFeatures.append(feature)

    return outputFeatures


def loadAnnotation(annotationFolder, identifier=None):
    annotationFiles = os.listdir(annotationFolder)
    annotationFiles = sorted([File
                              for File in annotationFiles
                              if File.endswith(".gbff")])

    if not annotationFiles:
        print("Annotation file not found! Check your annotation folder.")
        exit(1)

    annotationFile = annotationFiles[0]
    annotationFilePath = os.path.join(annotationFolder, annotationFile)
    annotationContent = list(SeqIO.parse(annotationFilePath, "genbank"))

    print(identifier)

    if identifier:
        wantedIdentifiers = [identifier, "chromosome_%s" % identifier]
    else:
        wantedIdentifiers = ["chromosome"]

    allIdentifiers = []

    noChromosomeIdentifierOnAnnotation = True
    for Scaffold in annotationContent:
        print(Scaffold.features[0].qualifiers)
        Qualifiers = Scaffold.features[0].qualifiers
        if 'chromosome' in Qualifiers.keys():
            noChromosomeIdentifierOnAnnotation = False
            ChromosomeName = Qualifiers['chromosome'][0]

            # Identifier identified;
            allIdentifiers.append(ChromosomeName)

            for wantedIdentifier in wantedIdentifiers:
                print(identifier)
                print(ChromosomeName)
                if wantedIdentifier.lower() == ChromosomeName.lower():
                    return Scaffold

    if noChromosomeIdentifierOnAnnotation and not identifier:
        print(">")
        return annotationContent[0]

    return allIdentifiers



