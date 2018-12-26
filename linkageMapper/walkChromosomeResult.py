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


class matrixViewer():
    def __init__(self, PWMData, drawPlotIndex):
        self.PWMData = PWMData
        self.drawPlot = drawPlotIndex

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
        self.figurecanvas.mpl_connect('button_press_event', self.nav_forward)
        self.index = 0

        vbox.pack_start(self.figurecanvas, True, True, 0)
        self.toolbar = NavigationToolbar(self.figurecanvas, win)
        self.toolbar.forward = self.nav_forward
        self.toolbar.back = self.nav_back
        self.toolbar.set_history_buttons()

        buttonBox = Gtk.HBox(True, 2)
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

    def nav_forward(self, d):
        self.change(self.index+1)
        self.index += 1

    def nav_back(self):
        pass

    def change(self, I):
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
        matrix, matrix_order, B = compute_serial_matrix(matrix, method="complete")
        labels = heatmapLabels[matrix_order]

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

    # plt.show()
    return fig


if __name__ == "__main__":
    last = None
    parser = OptionParser()

    parser.add_option("-i",
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

    # ITERATE PWM ANALYSIS DATA;
    viewer = matrixViewer(PWMData, plotPwmIndex)
    """
    for P in range(PWMData.shape[0]):

        if a == last:
            continue

        last = a
    """
