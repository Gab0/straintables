#!/bin/python

import argparse
import time
from Bio.Blast import NCBIWWW, NCBIXML
from straintables import OutputFile, Definitions


def ExtractAllPrimers(data):
    for i in range(data.content.shape[0]):
        for PT in Definitions.PrimerTypes:
            yield ("%s_%s" % (data.content.LocusName[i], PT), data.content[PT][i])


def Execute(options):

    PrimerData = OutputFile.MatchedRegions(options.WorkingDirectory)
    PrimerData.read()

    Query = [x[1] for x in (ExtractAllPrimers(PrimerData))]
    blast_result = NCBIWWW.qblast(program="blastn",
                                          database="nr", sequence=Query)


    T = 0
    while True:
        print(blast_result.read())
        time.sleep(1)
        T += 1
        contents = NCBIXML.parse(blast_result.read())
        if contents:
            print(contents)
            print("T=%i" % T)
            break


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", dest="WorkingDirectory")

    return parser.parse_args()


def main():
    options = parse_arguments()
    Execute(options)


if __name__ == "__main__":
    main()
