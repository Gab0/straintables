#!/bin/python

import pandas as pd
import os
import json

from . import Definitions


class OutputFile():
    def __init__(self, dirpath):
        self.dirpath = dirpath
        self.filepath = self.get_filepath()

    def get_filepath(self):
        return os.path.join(self.dirpath, self.filename)

    def check(self):
        return os.path.isfile(self.get_filepath())


class SimpleDataFrame():
    def add(self, data):
        self.content = pd.DataFrame(data, columns=self.columns)

    def write(self):
        self.content.to_csv(self.filepath, index=False)


class JsonFile():
    def write(self):
        with open(self.filepath, 'w') as f:
            json.dump(self.content, f, indent=2)

    def read(self):
        with open(self.filepath) as f:
            self.content = json.load(f)


class MatchedRegions(OutputFile, SimpleDataFrame):
    columns = [
        "LocusName",
        *Definitions.PrimerTypes,
        "RebootCount",
        "AlignmentHealth",
        "MeanLength",
        "StdLength",
        "Chromosome"
    ]
    filename = "MatchedRegions.csv"


class PrimerData(OutputFile, SimpleDataFrame):
    filename = "PrimerData.csv"


class AnalysisInformation(OutputFile, JsonFile):
    filename = "Information.json"
    fields = [
        "?"
    ]


class DockFailureReport(OutputFile, JsonFile):
    filename = "DockFailureReport.json"
