#!/bin/python

import unittest

from straintables.Database import StrainNames


class StrainNameExtractorTest(unittest.TestCase):
    def load_list(self, filepath):
        with open(filepath) as f:
            return [l.strip() for l in f.read().split("\n") if l]

    def test_extractor(self):
        Inputs = self.load_list("test/strain_headers/headers.txt")
        ExpectedOutputs = self.load_list("test/strain_headers/results.txt")

        for Input, Output in zip(Inputs, ExpectedOutputs):
            strain = StrainNames.fetchStrainName(Input)
            print("%s\n %s | %s\n" % (Input, strain, Output))
            assert strain == Output


if __name__ == "__main__":
    unittest.main()
