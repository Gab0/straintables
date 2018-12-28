#!/bin/python

import os
import pandas as pd
import numpy as np

from optparse import OptionParser

import matplotlib.pyplot as plt

import fastcluster
import scipy

import detectMutations

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvas
from matplotlib.backends.backend_gtk3 import (
    NavigationToolbar2GTK3 as NavigationToolbar)


# COMPLEX DISSIMILARITY MATRIX VIEWER GTK APPLICATION;
class matrixViewer():
    def __init__(self, PWMData, drawPlotIndex, allowedIndexes):
        self.PWMData = PWMData
        self.drawPlot = drawPlotIndex
        self.allowedIndexes = allowedIndexes

        win = Gtk.Window()
        win.connect("destroy", lambda x: Gtk.main_quit())
        win.set_default_size(400, 300)
        win.set_title("Embedding in GTK")

        vbox = Gtk.VBox()
        win.add(vbox)

        # INITIALIZE PLOT FIGURE;
        self.figure = plt.figure()

        self.drawPlot(self.figure, self.PWMData, 0)

        self.figurecanvas = FigureCanvas(self.figure)  # a Gtk.DrawingArea
        self.index = 0

        vbox.pack_start(self.figurecanvas, True, True, 0)
        self.toolbar = NavigationToolbar(self.figurecanvas, win)
        self.toolbar.forward = self.nav_forward
        self.toolbar.back = self.nav_back
        self.toolbar.set_history_buttons()

        buttonBox = Gtk.HBox(homogeneous=True, spacing=2)
        self.btn_back = Gtk.Button(label="<")
        self.btn_back.connect("clicked", self.nav_back)
        self.btn_next = Gtk.Button(label=">")
        self.btn_next.connect("clicked", self.nav_forward)

        buttonBox.add(self.btn_back)
        buttonBox.add(self.btn_next)

        vbox.pack_start(buttonBox, False, False, 0)
        vbox.pack_start(self.toolbar, False, False, 0)

        win.show_all()
        Gtk.main()

    def cycleIndexes(self, amt):
        self.index += amt

        while self.index not in self.allowedIndexes:
            self.index += amt
            if self.index > max(self.allowedIndexes):
                self.index = min(self.allowedIndexes)

    def nav_forward(self, d):
        self.cycleIndexes(1)
        self.changeView(self.index)

    def nav_back(self, d):
        self.cycleIndexes(-1)
        self.changeView(self.index)

    def changeView(self, I):
        self.figure.clf()
        self.drawPlot(self.figure, self.PWMData, I)
        self.figurecanvas.draw()
        self.figurecanvas.flush_events()


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


def reorderMatrix(original_matrix, res_order):
    N = len(original_matrix)
    reordered_matrix = np.zeros((N, N))

    a, b = np.triu_indices(N, k=1)
    reordered_matrix[a, b] = original_matrix[[res_order[i] for i in a], [res_order[j] for j in b]]
    reordered_matrix[b, a] = reordered_matrix[a, b]
    return reordered_matrix


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

    seriated_dist = reorderMatrix(dist_mat, res_order)

    return seriated_dist, res_order, res_linkage


def createSubplot(fig, position, name, matrix, labels):
    new_ax = fig.add_subplot(position)

    detectMutations.heatmapToAxis(matrix, new_ax, labels=labels)

    new_ax.set_xlabel(name)


def plotPwmIndex(fig, PWMData, I):
    d = PWMData.iloc[I]
    a = d["Unnamed: 0"]
    b = d["Unnamed: 1"]

    # walk loci by loci mode.

    # EXTRACR LOCUS NAMES;
    a_name, b_name = fixArrayFilename(a), fixArrayFilename(b)

    try:
        data = [
            PrimerData[PrimerData.Locus == name.replace("LOCI_", "")].iloc[0]
            for name in [a_name, b_name]
        ]
    except IndexError:
        print("Failure on %s" % a_name)

    # LOAD MATRIX DATA;
    ma = np.load(buildArrayPath(a))
    mb = np.load(buildArrayPath(b))

    # REORDERED MATRIXES;
    ordered_ma, matrix_order, B = compute_serial_matrix(ma, method="complete")
    ordered_mb = reorderMatrix(mb, matrix_order)
    orderedLabels = heatmapLabels[matrix_order]

    # plot;
    createSubplot(fig, 331, a_name, ordered_ma, orderedLabels)
    createSubplot(fig, 333, b_name, ordered_mb, orderedLabels)

    # ORIGINAL MATRIXES;
    # plot;
    createSubplot(fig, 337, a_name, ma, heatmapLabels)
    createSubplot(fig, 339, b_name, mb, heatmapLabels)

    # BUILD SHOWN INFO;
    distance = abs(data[0].PositionStart - data[1].PositionStart)
    Title = [
        "Distance = {:,} bp".format(distance),
        "%s vs %s" % (a_name, b_name),
        "Mantel=%.4f     p=%.4f" % (d["mantel"], d["mantel_p"]),
        "DIFF=%i" % d["matrix_ranking_diff"],
        " "
    ]

    Title = "\n".join(Title)

    # ADDITIONAL INFORMATION FIGURE;
    ax_t = fig.add_subplot(335)

    ax_t.text(-0.2,
              0.6,
              s=Title,
              clip_on=False
    )

    ax_t.axis("off")

    plt.title("")

    return fig


if __name__ == "__main__":
    last = None
    parser = OptionParser()

    parser.add_option("-d",
                      dest="inputDirectory")

    options, args = parser.parse_args()

    PrimerFilePath = os.path.join(options.inputDirectory, "PrimerData.csv")
    PWMFilePath = os.path.join(options.inputDirectory, "PWMAnalysis.csv")

    PrimerData = pd.read_csv(PrimerFilePath)
    PWMData = pd.read_csv(PWMFilePath)

    # FETCH ORIGINAL HEATMAP GENOME LABELS;
    heatmapLabelsFilePath = os.path.join(
        options.inputDirectory,
        "heatmap_labels.npy"
    )
    heatmapLabels = np.load(heatmapLabelsFilePath)

    # FETCH VIEWABLE DATA INDEXES;
    allowedIndexes = []
    for I in range(PWMData.shape[0]):
        d = PWMData.iloc[I]
        a = d["Unnamed: 0"]
        if a == last:
            continue
        allowedIndexes.append(I)
        last = a

    # SHOW DATA;
    viewer = matrixViewer(PWMData, plotPwmIndex, allowedIndexes)
