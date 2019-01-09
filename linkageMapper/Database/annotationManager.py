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
    annotationFiles = sorted([file for file in annotationFiles if file.endswith(".gbff")])

    if not annotationFiles:
        print("Annotation file not found! Check your annotation folder.")
        exit(1)

    annotationFile = annotationFiles[0]
    annotationFilePath = os.path.join(annotationFolder, annotationFile)
    annotationContent = SeqIO.parse(annotationFilePath, "genbank")

    for Scaffold in annotationContent:
        ChromosomeName = Scaffold.features[0].qualifiers['chromosome'][0]
        if identifier.lower() == ChromosomeName.lower():
            return Scaffold

    return None

