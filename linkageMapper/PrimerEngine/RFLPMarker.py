#!/bin/python

import os
import pandas as pd


class RFLPReference():
    def getDatabase(self, filename):
        filepath = os.path.join(self.inputDirectory, filename)
        return pd.read_csv(filepath)

    def __init__(self, inputDirectory):
        self.genotypeData = self.getDatabase("genomes_haplogroups.csv")
        self.rflpGenotypes = self.getDatabase("genotypes.csv")

    def getGenotypeNumber(self, name):
        found = self.genotypeData[self.genotypeData.Genome == name]
        GenotypeNumber = found.iloc[0].ToxoDB

        return GenotypeNumber

    def getRFLPLocus(self, genotypeNumber, referenceLocus):

        found = self.rflpGenotypes[
            self.rflpGenotypes.Genotype == genotypeNumber]
        RFLPLocus = found.iloc[0][referenceLocus]
        return RFLPLocus


