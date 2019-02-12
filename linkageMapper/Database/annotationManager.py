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

    allIdentifiers = []
    for Scaffold in annotationContent:
        Qualifiers = Scaffold.features[0].qualifiers
        if 'chromosome' in Qualifiers.keys():
            ChromosomeName = Qualifiers['chromosome'][0]
            allIdentifiers.append(ChromosomeName)
            if identifier:
                if identifier.lower() == ChromosomeName.lower():
                    return Scaffold

    if not identifier:
        return allIdentifiers



