#!/bin/python

import pandas as pd
import os

from . import Definitions


class SimpleDataFrame():
    def __init__(self, data):
        self.content = pd.DataFrame(data, columns=self.columns)

    def write(self, dirpath):
        filepath = os.path.join(dirpath, self.filename)
        self.content.to_csv(filepath, index=False)


class MatchedPrimers(SimpleDataFrame):
    columns = [
        "LocusName",
        *Definitions.PrimerTypes,
        "RebootCount",
        "AlignmentHealth",
        "MeanLength",
        "StdLength"
    ]
    filename = "MatchedRegions.csv"
