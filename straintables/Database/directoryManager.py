#!/bin/python

import os


def createDirectoryPath(DirectoryPath):
    if DirectoryPath and not os.path.isdir(DirectoryPath):
        Path = [
            step for step in DirectoryPath.split(os.sep)
            if step
        ]

        for d, Directory in enumerate(Path):
            subDirectoryPath = os.path.join(*Path[:d+1])
            if not os.path.isdir(subDirectoryPath):
                os.mkdir(subDirectoryPath)
                print("Creating directory %s." % subDirectoryPath)
