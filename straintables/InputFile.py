#!/bin/python

import pandas as pd
import numpy as np


# LOAD USER DEFINED PRIMER DATA, with or without header;
def loadPrimerList(filePath):

    lociPrimerList = pd.read_csv(filePath)

    expectedColumns = ["LocusName", "ForwardPrimer", "ReversePrimer"]

    fileColumns = list(lociPrimerList.columns)

    for i in range(len(fileColumns)):
        if "Unnamed" in fileColumns[i]:
            fileColumns[i] = np.nan

    if fileColumns != expectedColumns:
        newFirstRowData = dict([(expected, fileColumns[e])
                                for e, expected
                                in enumerate(expectedColumns)])

        newFirstRow = pd.DataFrame([newFirstRowData],
                                   columns=expectedColumns)

        lociPrimerList.columns = expectedColumns

        lociPrimerList = pd.concat([newFirstRow, lociPrimerList],
                                   axis=0,
                                   ignore_index=True).reset_index(drop=True)

    return lociPrimerList
