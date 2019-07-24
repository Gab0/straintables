#!/bin/python

import pandas as pd

from . import Definitions


class MatchedPrimers():
    columns = [
        "LocusName",
        *Definitions.PrimerTypes,
        "RebootCount",
        "AlignmentHealth",
        "MeanLength",
        "StdLength"
    ]

    def __init__(self, data):
        self.content = pd.DataFrame(data, columns=self.columns)

    def write(self, filepath):
        self.content.to_csv(filepath, index=False)
