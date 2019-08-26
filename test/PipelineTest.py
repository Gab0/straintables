#!/bin/python

import unittest
import os
import shutil

from straintables.Executable import Pipeline, fetchDataNCBI, MatrixViewer


class PipelineTest(unittest.TestCase):
    WorkingDirectory = "straintables_test"
    PrimerFile = "t.csv"
    Regions = ["moeZ", "gmk", "pntB"]
    PrimerFilePath = os.path.join(WorkingDirectory, PrimerFile)

    def step1_prepare_directory(self):
        print("Initializing test directory at %s..." % self.WorkingDirectory)
        if os.path.isdir(self.WorkingDirectory):
            shutil.rmtree(self.WorkingDirectory)

        os.mkdir(self.WorkingDirectory)

        with open(self.PrimerFilePath, 'w') as f:
            f.write("\n".join(self.Regions))

    def step2_download(self):
        download_options = fetchDataNCBI.parse_arguments()
        download_options.queryOrganism = "Mycobacterium leprae"
        download_options.WorkingDirectory = self.WorkingDirectory
        fetchDataNCBI.Execute(download_options)

    def step3_run_pipeline(self):

        pipeline_options = Pipeline.parse_arguments()
        pipeline_options.WorkingDirectory = self.WorkingDirectory
        pipeline_options.PrimerFile = self.PrimerFilePath
        pipeline_options.SourceDataDirectory = self.WorkingDirectory
        Pipeline.Execute(pipeline_options)

    def step4_run_viewer(self):

        viewer_options = MatrixViewer.parse_arguments()
        viewer_options.inputDirectory = self.WorkingDirectory
        viewer_options.inputDir = self.WorkingDirectory

        MatrixViewer.Execute(viewer_options)

    def _steps(self):
        for name in dir(self):
            if name.startswith("step"):
                yield name, getattr(self, name)

    def test_steps(self):
        for name, step in self._steps():
            try:
                step()
            except Exception as e:
                self.fail("{} failed ({}: {})".format(name, str(type(e)), e))


if __name__ == "__main__":
    unittest.main()
