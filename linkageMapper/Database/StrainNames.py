#!/bin/python

import re


def fetchStrainName(genomeDescriptor, organismName=None):

    strainName = None

    # RENAMING CRITERIA 1:
    def renamingCriteria1():
        queries = [
            "strain ([\w]+ [\d]+)"
            "strain ([A-Z\-0-9]+)",
            "[\w]+ [\d]+"
            "[\w\-\d]+",
            #"%s strain ([\w]+ [\d]+)" % organismName,
            #"%s strain ([\w\-\d]+)" % organismName,

        ]
        for q, Query in enumerate(queries):
            d = re.findall(Query,
                           genomeDescriptor
)
            if d:
                print(q)
                return d[0]

    # RENAMING CRITERIA 2:
    def renamingCriteria2():
        e = re.findall("%s ([\w\-\d]+)" % (organismName),
                       genomeDescriptor,
                       flags=re.IGNORECASE)
        if e:
            return e[0]

    # RENAMING CRITERIA 3:
    def renamingCriteria3():
        words = genomeDescriptor.split(" ")[1:]
        uppercaseWords = [w.isupper() for w in words]
        if sum(uppercaseWords) == 1:
            return words[uppercaseWords.index(True)]

    allCriteria = [
        renamingCriteria1,
        #renamingCriteria2,
        renamingCriteria3
    ]

    genomeDescriptor = genomeDescriptor.replace("isolate ", "")
    for renamingCriteria in allCriteria:
        strainName = renamingCriteria()
        if strainName is not None:
            print(renamingCriteria.__name__)
            return strainName.replace(" ", "_")

    return None
