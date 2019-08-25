#!/bin/python

import pandas as pd
import os
import json

from . import Definitions


# -- BASE CLASS FOR ALL TYPES OF OUTPUT FILES;
class OutputFile():
    def __init__(self, dirpath):
        self.dirpath = dirpath
        self.filepath = self.get_filepath()

    def get_filepath(self):
        return os.path.join(self.dirpath, self.filename)

    def check(self):
        return os.path.isfile(self.get_filepath())


# -- TYPES OF OUTPUT FILES;
class SimpleDataFrame(OutputFile):
    def add(self, data):
        self.content = pd.DataFrame(data, columns=self.columns)

    def write(self):
        self.content.to_csv(self.filepath, index=False)

    def read(self):
        self.content = pd.read_csv(self.filepath)

class JsonFile(OutputFile):
    content = {}

    def write(self):
        with open(self.filepath, 'w') as f:
            json.dump(self.content, f, indent=2)

    def read(self):
        with open(self.filepath) as f:
            self.content = json.load(f)


# -- OUTPUT FILE FLAVORS;
class MatchedRegions(SimpleDataFrame):
    columns = [
        "LocusName",
        *Definitions.PrimerTypes,
        "RebootCount",
        "AlignmentHealth",
        "MeanLength",
        "StdLength",
        "Chromosome",
        "StartPosition"
    ]
    filename = "MatchedRegions.csv"


class PrimerData(SimpleDataFrame):
    filename = "PrimerData.csv"


class AnalysisInformation(JsonFile):
    filename = "Information.json"
    fields = [
        "?"
    ]


class DockFailureReport(JsonFile):
    filename = "DockFailureReport.json"
