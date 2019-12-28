#!/bin/python

import argparse
import time
from Bio.Blast import NCBIWWW, NCBIXML
from straintables import OutputFile, Definitions


def Execute(options):

    data = OutputFile.MatchedRegions(options.WorkingDirectory)
    data.read()

    def ExtractFromData(data):
        for i in range(data.content.shape[0]):
            for PT in Definitions.PrimerTypes:
                yield data.content[PT][i]

    Query = list(ExtractFromData(data))
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
