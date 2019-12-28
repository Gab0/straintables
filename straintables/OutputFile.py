#!/bin/python

import pandas as pd
import numpy as np
import os
import json

from . import Definitions


# -- BASE CLASS FOR ALL TYPES OF OUTPUT FILES;
class OutputFile():
    def __init__(self, dirpath):
        self.dirpath = dirpath
        self.filepath = self.get_filepath()
        # if self.check():
        #    self.read()

    def get_filepath(self):
        return os.path.join(self.dirpath, self.filename)

    def check(self):
        return os.path.isfile(self.get_filepath())


# -- DATA FRAME OUTPUT FILE;
class SimpleDataFrame(OutputFile):
    csv_index = False

    def add(self, data):
        self.content = pd.DataFrame(data, columns=self.columns)

    def write(self):
        self.content.to_csv(self.filepath, index=self.csv_index)

    def read(self):
        self.content = pd.read_csv(self.filepath)


# -- JSON OUTPUT FILE;
class JsonFile(OutputFile):
    content = {}

    def write(self):
        with open(self.filepath, 'w') as f:
            json.dump(self.content, f, indent=2)

    def read(self):
        with open(self.filepath) as f:
            self.content = json.load(f)


# -- OUTPUT FILE INSTANTIABLE CLASSES;
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
    columns = ["Locus",
               "Sequence",
               "PositionStart",
               "PositionEnd"
               ]
    filename = "PrimerData.csv"


class AlignedRegions(SimpleDataFrame):
    filename = "AlignedRegions.csv"


class AnalysisInformation(JsonFile):
    filename = "Information.json"
    fields = [
        "? - TBD"
    ]


class DockFailureReport(JsonFile):
    filename = "DockFailureReport.json"


class PWMAnalysis(SimpleDataFrame):
    csv_index = True
    filename = "PWMAnalysis.csv"

# WIP
class ValueMatrix(OutputFile):
    def save(self):
        np.save()
