#!/bin/python
import random
import pandas as pd

from optparse import OptionParser

import annotationManager

parser = OptionParser()

parser.add_option("-i",
                  dest="inputAnnotation")

parser.add_option("-o",
                  dest="outputFile")


options, args = parser.parse_args()


allFeatures = annotationManager.loadFeatures(options.inputAnnotation)

chosenFeatures = []
for feature in allFeatures:
    COND1 = "gene" in feature.qualifiers.keys()
    COND2 = random.random() < 0.03
    NAME = None
    if COND1:
        NAME = feature.qualifiers["gene"][0]
    elif COND2:
        NAME = feature.qualifiers["locus_tag"][0]

    if NAME:
        chosenFeatures.append(NAME)


data = [{"LocusName": feature} for feature in chosenFeatures]
data = pd.DataFrame(data, columns=[
    "LocusName",
    "ForwardPrimer",
    "ReversePrimer"
]
)
data.to_csv(options.outputFile, index=False)
