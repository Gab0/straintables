#!/bin/python 

import ete3
import sys

if len(sys.argv) < 2:
    exit("No file specified.")

filePath = sys.argv[1]


t = ete3.Tree(filePath)
ts = ete3.TreeStyle()
ts.show_branch_length = True
t.render(filePath + ".pdf", h=500, tree_style=ts)
