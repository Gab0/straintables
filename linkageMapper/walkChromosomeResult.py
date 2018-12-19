#!/bin/python

import os
import pandas as pd
import numpy as np

from optparse import OptionParser

import matplotlib.pyplot as plt

import detectMutations

parser = OptionParser()

parser.add_option("-i",
                  dest="inputDirectory")

options, args = parser.parse_args()


PrimerFilePath = os.path.join(options.inputDirectory, "PrimerData.csv")
PWMFilePath = os.path.join(options.inputDirectory, "PWMAnalysis.csv")

PrimerData = pd.read_csv(PrimerFilePath)
PWMData = pd.read_csv(PWMFilePath)


def buildArrayPath(f):
    return os.path.join(options.inputDirectory, f)


def fixArrayFilename(f):
    return f.split('.')[0]


last = None
heatmapLabels = np.load(os.path.join(options.inputDirectory, "heatmap_labels.npy"))
for P in range(PWMData.shape[0]):

    d = PWMData.iloc[P]
    a = d["Unnamed: 0"]
    b = d["Unnamed: 1"]

    # walk loci by loci mode.
    if a == last:
        continue

    # EXTRACR LOCUS NAMES;
    a_name, b_name = fixArrayFilename(a), fixArrayFilename(b)

    # LOAD MATRIX DATA;
    ma = np.load(buildArrayPath(a))
    mb = np.load(buildArrayPath(b))

    # labels = buildArrayPath("labels.")
    fig = plt.figure()

    ax_a = fig.add_subplot(212)
    detectMutations.heatmapToAxis(ma, ax_a, labels=heatmapLabels)
    ax_a.set_title(a_name, loc='left', pad=-240)

    ax_b = fig.add_subplot(221)
    detectMutations.heatmapToAxis(mb, ax_b, labels=heatmapLabels)
    ax_b.set_title(b_name, loc='left', pad=-240)

    try:
        data = [
            PrimerData[PrimerData.Locus == name.replace("LOCI_", "")].iloc[0]
            for name in [a_name, b_name]
        ]
    except IndexError:
        print("Failure on %s" % a_name)
        continue

    Title = [
        "Distance = %ibp" % (abs(data[0].PositionStart - data[1].PositionStart)),
        "%s vs %s" % (a_name, b_name),
        "Mantel=%.4f     p=%.4f" % (d["mantel"], d["mantel_p"]),
        "DIFF=%i" % d["matrix_ranking_diff"],
        " "
    ]

    Title = "\n".join(Title)

    ax_t = fig.add_subplot(222)

    ax_t.text(-0.2,
              0.6,
              s=Title,
              clip_on=False
    )
    ax_t.axis("off")

    plt.title("")
    # plt.tight_layout()

    plt.show()

    last = a


"""
for P in range(PrimerData.shape[0]):
    if not P % 2:
        Primer = PrimerData.iloc[P]
        matrixFile = os.path.join(options.inputDirectory,
                                  "LOCI_%s.aln.pdf" % Primer["Locus"])

        subprocess.run(["okular", matrixFile])
"""
