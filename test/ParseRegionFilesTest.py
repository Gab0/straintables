#!/bin/python

import unittest
import os
from straintables import InputFile


class ParseRegionFilesTest(unittest.TestCase):
    test_files_dir = "test/region_files"
    def test_inputs(self):

        InputFiles = [
            os.path.join(self.test_files_dir, f)
            for f in os.listdir(self.test_files_dir)
        ]

        for Input in InputFiles:
            data = InputFile.loadPrimerList(Input)
            print(data)
            assert list(data["LocusName"]) == ["R1", "R2", "r3"]


if __name__ == "__main__":
    unittest.main()
