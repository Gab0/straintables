#!/bin/python

import pandas as pd
import numpy as np


# LOAD USER DEFINED PRIMER DATA, with or without header;
def loadPrimerList(filePath):

    lociPrimerList = pd.read_csv(filePath)

    expectedColumns = ["LocusName", "ForwardPrimer", "ReversePrimer"]

    fileColumns = list(lociPrimerList.columns)

    while len(fileColumns) < len(expectedColumns):
        fileColumns.append(np.nan)

    if fileColumns != expectedColumns:
        lociPrimerList = pd.read_csv(filePath, header=None)

        # -- Make sure data has the correct number of columns;
        while lociPrimerList.shape[1] < len(expectedColumns):
            lociPrimerList[str(lociPrimerList.shape[1])] =\
                np.full(lociPrimerList.shape[0], np.nan)

        if lociPrimerList.shape[1] > len(expectedColumns):
            lociPrimerList = lociPrimerList.loc[:, :len(expectedColumns)-1]

        lociPrimerList.columns = expectedColumns

    return lociPrimerList
