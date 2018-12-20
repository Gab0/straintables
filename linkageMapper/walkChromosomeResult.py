#!/bin/python

import os
import pandas as pd
import numpy as np

from optparse import OptionParser

import matplotlib.pyplot as plt

import fastcluster
import scipy

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


def seriation(Z, N, cur_index):
    '''
        input:
            - Z is a hierarchical tree (dendrogram)
            - N is the number of points given to the clustering process
            - cur_index is the position in the tree for the recursive traversal
        output:
            - order implied by the hierarchical tree Z

        seriation computes the order implied by a hierarchical tree (dendrogram)
    '''
    if cur_index < N:
        return [cur_index]
    else:
        left = int(Z[cur_index - N,0])
        right = int(Z[cur_index - N, 1])
        return (seriation(Z, N, left) + seriation(Z, N, right))


def compute_serial_matrix(dist_mat, method="ward"):
    '''
        input:
            - dist_mat is a distance matrix
            - method = ["ward","single","average","complete"]
        output:
            - seriated_dist is the input dist_mat,
              but with re-ordered rows and columns
              according to the seriation, i.e. the
              order implied by the hierarchical tree
            - res_order is the order implied by
              the hierarhical tree
            - res_linkage is the hierarhical tree (dendrogram)

        compute_serial_matrix transforms a distance matrix into 
        a sorted distance matrix according to the order implied 
        by the hierarchical tree (dendrogram)
    '''
    N = len(dist_mat)
    flat_dist_mat = scipy.spatial.distance.squareform(dist_mat)
    res_linkage = fastcluster.linkage(flat_dist_mat, method=method, preserve_input=True)
    res_order = seriation(res_linkage, N, N + N-2)
    seriated_dist = np.zeros((N, N))
    a, b = np.triu_indices(N, k=1)
    seriated_dist[a, b] = dist_mat[ [res_order[i] for i in a], [res_order[j] for j in b]]
    seriated_dist[b, a] = seriated_dist[a, b]

    return seriated_dist, res_order, res_linkage


def createSubplot(fig, position, name, matrix, labels, Reorder):
    new_ax = fig.add_subplot(position)
    if Reorder:
        # REORDER MATRIX
        matrix, matrix_order, B = compute_serial_matrix(1-matrix, method="complete")
        labels = heatmapLabels[matrix_order]

    detectMutations.heatmapToAxis(matrix, new_ax, labels=labels)

    new_ax.set_xlabel(name)


if __name__ == "__main__":
    last = None

    # FETCH ORIGINAL HEATMAP GENOME LABELS;
    heatmapLabelsFilePath = os.path.join(
        options.inputDirectory,
        "heatmap_labels.npy"
    )
    heatmapLabels = np.load(heatmapLabelsFilePath)

    # ITERATE PWM ANALYSIS DATA;
    for P in range(PWMData.shape[0]):

        d = PWMData.iloc[P]
        a = d["Unnamed: 0"]
        b = d["Unnamed: 1"]

        # walk loci by loci mode.
        if a == last:
            continue

        # EXTRACR LOCUS NAMES;
        a_name, b_name = fixArrayFilename(a), fixArrayFilename(b)

        try:
            data = [
                PrimerData[PrimerData.Locus == name.replace("LOCI_", "")].iloc[0]
                for name in [a_name, b_name]
            ]
        except IndexError:
            print("Failure on %s" % a_name)
            continue

        # LOAD MATRIX DATA;
        ma = np.load(buildArrayPath(a))
        mb = np.load(buildArrayPath(b))

        # INITIALIZE PLOT FIGURE;
        fig = plt.figure()

        # ORIGINAL MATRIXES;
        createSubplot(fig, 331, a_name, ma, heatmapLabels, True)
        createSubplot(fig, 333, b_name, mb, heatmapLabels, True)

        # REORDERED MATRIXES;
        createSubplot(fig, 337, a_name, ma, heatmapLabels, False)
        createSubplot(fig, 339, b_name, mb, heatmapLabels, False)

        # BUILD SHOWN INFO;
        Title = [
            "Distance = %ibp" % (abs(data[0].PositionStart - data[1].PositionStart)),
            "%s vs %s" % (a_name, b_name),
            "Mantel=%.4f     p=%.4f" % (d["mantel"], d["mantel_p"]),
            "DIFF=%i" % d["matrix_ranking_diff"],
            " "
        ]

        Title = "\n".join(Title)

        ax_t = fig.add_subplot(335)

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
