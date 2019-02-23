#!/bin/python
import random
import pandas as pd

from optparse import OptionParser

from Database import annotationManager


def Execute(options):
    selectedScaffold = annotationManager.loadAnnotation(
        options.inputAnnotationFolder,
        identifier=options.inputAnnotationName
    )

    """
    if not options.inputAnnotationName:
        print("No chromosome selected. Available Chromosomes:")
        availableChromosomes = [c for c in selectedScaffold if "Unknown" not in c]
        print("\n".join(availableChromosomes))
        return
    """

    if selectedScaffold is None or type(selectedScaffold) == list:

        print("Chromosome %s not found." % options.inputAnnotationName)
        print(selectedScaffold)
        return

    print("Found feature source scaffold, please review: \n\n%s" % selectedScaffold)
    chosenFeatures = []
    for feature in selectedScaffold.features:
        COND1 = "gene" in feature.qualifiers.keys()
        COND2 = "locus_tag" in feature.qualifiers.keys() and random.random() < options.locusProbability

        NAME = None
        if COND1:
            NAME = feature.qualifiers["gene"][0]
        elif COND2:
            NAME = feature.qualifiers["locus_tag"][0]

        if NAME:
            chosenFeatures.append(NAME)

    # IT REPEATS A LOT OF LOCI, BECAUSE THE ANNOTATION HAS MANY ENTRIES WITH SAME NAME.
    # MAYBE MAKING A SET IS SUFFICIENT?
    chosenFeatures = list(set(chosenFeatures))

    data = [{"LocusName": feature} for feature in chosenFeatures]
    data = pd.DataFrame(data, columns=[
        "LocusName",
        "ForwardPrimer",
        "ReversePrimer"
    ]
    )

    data.to_csv(options.outputFile, index=False)

    print("\n%s written with %i primers." % (options.outputFile, data.shape[0]))


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-d",
                      dest="inputAnnotationFolder",
                      default='annotations')

    parser.add_option("-c", "--chromosome",
                      dest="inputAnnotationName",
                      help="Chromosome identifier. E.g 'X' or 'II' or 'Ia'")

    parser.add_option("-o",
                      dest="outputFile")

    parser.add_option("-p",
                      dest="locusProbability",
                      type="float",
                      default=0.1)

    parser.add_option('-l',
                      dest="ListChromosomeNames",
                      action="store_true",
                      help="Print all chromosome names and exit.")

    options, args = parser.parse_args()

    print("\n\n")
    Execute(options)
