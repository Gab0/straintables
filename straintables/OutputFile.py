#!/bin/python

import pandas as pd
import os
import json

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


class AnalysisInformation():
    filename = "Information.json"
    fields = [
        "?"
    ]

    def write(self, dirpath):
        filepath = os.path.join(dirpath, self.filename)
        with open(filepath, 'w') as f:
            json.dump(self.content, f, indent=2)
