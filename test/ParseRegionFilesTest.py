#!/bin/python

import unittest
import os
from straintables import InputFile


class ParseRegionFilesTest(unittest.TestCase):
    def TestInputs(self):

        InputFiles = [
            os.path.join("test/region_files", f)
            for f in os.listdir("region_files")
        ]

        for Input in InputFiles:
            assert InputFile.loadPrimerList(Input) == ["GA", "Ga", "ga"]


if __name__ == "__main__":
    unittest.main()
