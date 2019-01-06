#!/bin/python
import re
import os
import pandas as pd
import numpy as np

import subprocess

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

        vbox.pack_start(self.figurecanvas, True, True, 0)

        # INITIALIZE STATE;
        self.labelColorsOn = 1
        self.index = 0


        # SHOW LOCUS NAVIGATION TOOLBAR;
        buttonBox = Gtk.HBox(homogeneous=True, spacing=2)
        self.btn_back = Gtk.Button(label="<")
        self.btn_back.connect("clicked", self.nav_back)
        self.btn_next = Gtk.Button(label=">")
        self.btn_next.connect("clicked", self.nav_forward)

        buttonBox.add(self.btn_back)
        buttonBox.add(self.btn_next)

        vbox.pack_start(buttonBox, False, False, 0)

        # MODIFY MATPLOTLIB TOOLBAR;
        self.toolbar = NavigationToolbar(self.figurecanvas, win)
        self.toolbar.forward = self.nav_forward
        self.toolbar.back = self.nav_back
        self.toolbar.set_history_buttons()

        # SET BOTTOM TOOLBAR, WHICH INCLUDE MATPLOTLIB BAR;
        panelBox = Gtk.HBox(homogeneous=False, spacing=2)
        toggleColor = Gtk.ToggleButton(label=None, image=Gtk.Image(stock=Gtk.STOCK_COLOR_PICKER))
        toggleColor.set_tooltip_text("Show Matrix Label Colors")
        toggleColor.set_active(True)
        toggleColor.connect("clicked", self.toggleColor)

        infoText = Gtk.Label.new(options.inputDirectory)

        panelBox.pack_start(self.toolbar, expand=False, fill=False, padding=0)
        panelBox.pack_start(toggleColor, False, False, 0)

        panelBox.add(infoText)

        vbox.pack_start(panelBox, False, False, 0)

        cid = self.figurecanvas.mpl_connect('button_press_event', self.onclickCanvas)

        win.show_all()

        # HIDE THIS ANNOYING THING.
        self.toolbar.message.hide()

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

    def toggleColor(self, d):
        self.labelColorsOn = 1 - self.labelColorsOn
        self.changeView(self.index)

    def changeView(self, I):
        self.figure.clf()
        self.drawPlot(self.figure, self.PWMData, I, showLabelColors=self.labelColorsOn)
        self.figurecanvas.draw()
        self.figurecanvas.flush_events()

    def onclickCanvas(self, event):
        # print(event)
        if event.inaxes:
            Data = self.PWMData.iloc[self.index]

            bbox = self.figure.get_window_extent().transformed(
                self.figure.dpi_scale_trans.inverted())
            width, height = bbox.width*self.figure.dpi, bbox.height*self.figure.dpi

            if event.x < width // 2:
                KeyName = "Unnamed: 0"
            else:
                KeyName = "Unnamed: 1"

            LocusName = Data[KeyName]
            alignmentFilePath = os.path.join(options.inputDirectory,
                                             LocusName.replace(".npy", ""))
            command = ["aliview", alignmentFilePath]
            print(command)

            subprocess.run(command)


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
    reordered_matrix[a, b] = original_matrix[[res_order[i] for i in a],
                                             [res_order[j] for j in b]]
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

    return new_ax


def singleLocusStatus(axis, locus_name):

    # FETCH HEALTH SCORE FOR LOCUS;
    locus_identifier = locus_name.replace("LOCI_", "")
    Health = MatchData[MatchData.LocusName == locus_identifier]
    if not Health.empty:
        Health = Health.iloc[0]["AlignmentHealth"]

    # DECLARE DISPLAY COLORS;
    colorRanges = {
        "red": (0, 50),
        "orange": (50, 70),
        "green": (70, 100)
    }

    # SELECT DISPLAY COLORS;
    color = "black"
    for anycolor in colorRanges.keys():
        v = colorRanges[anycolor]
        if v[0] <= Health <= v[1]:
            color = anycolor

    # PRINT ADJACENT TEXT;
    axis.text(-0.2,
              0.6,
              s="Amplicon Health:",
              clip_on=False,
              fontsize=12)

    # PRINT COLORED HEALTH VALUE TEXT;
    axis.text(0.4,
              0.6,
              s="%.2f%%" % Health,
              clip_on=False,
              color=color,
              fontsize=15)

    # DISABLE AXIS XY RULERS;
    axis.axis("off")


def parseMeshcluster(clusterFilePath):
    clusterData = open(clusterFilePath).read().split("\n")
    clusterOutputData = {}
    for line in clusterData:
        if ">Cluster" in line:
            key = int(re.findall(">Cluster (\d+)", line)[0])
        if "nt," in line:
            if key not in clusterOutputData.keys():
                clusterOutputData[key] = []
            Individual = re.findall(">([\w\d]+)...", line)[0]
            clusterOutputData[key].append(Individual)
    # print(clusterOutputData)
    return clusterOutputData


def plotPwmIndex(fig, PWMData, I, showLabelColors=True):
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
    r_axis1 = createSubplot(fig, 331, a_name, ordered_ma, orderedLabels)
    r_axis2 = createSubplot(fig, 333, b_name, ordered_mb, orderedLabels)


    # ORIGINAL MATRIXES;
    # plot;
    o_axis1 = createSubplot(fig, 337, a_name, ma, heatmapLabels)
    o_axis2 = createSubplot(fig, 339, b_name, mb, heatmapLabels)


    # COLORIZE MATRIX LABELS BY MESHCLUSTER;
    if showLabelColors:
        plotFilenames = [
            (r_axis1, a_name),
            (r_axis2, b_name),
            (o_axis1, a_name),
            (o_axis2, b_name)
        ]

        colorMap = plt.get_cmap("Set1")

        for (Axis, LocusName) in plotFilenames:
            clusterFilePath = buildArrayPath(LocusName) + ".clst"
            if os.path.isfile(clusterFilePath):
                clusterOutputData = parseMeshcluster(clusterFilePath)
                NB_Groups = len(clusterOutputData.keys())
                GroupColors = [colorMap(x / NB_Groups)
                               for x in range(NB_Groups)]

                # COLORIZE LABELS;
                for (labelx, labely) in zip(Axis.get_xticklabels(),
                                            Axis.get_yticklabels()):
                    text = labelx.get_text()
                    for key in clusterOutputData.keys():
                        if text in clusterOutputData[key]:
                            labelx.set_color(GroupColors[key])
                            labely.set_color(GroupColors[key])
                            break

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

    # ALIGNMENT HEALTH INFORMATION FIGURE;
    if "AlignmentHealth" in MatchData.keys():
        ax_ha = fig.add_subplot(334)
        ax_hb = fig.add_subplot(336)

        singleLocusStatus(ax_ha, a_name)
        singleLocusStatus(ax_hb, b_name)

        # Additional info on secondary axis
        RecombinationMessage = "True" if d["recombination"] else "False"
        Message = "Recombination? %s" % RecombinationMessage
        ax_hb.text(0.8, 1, s=Message)

    # RECOMBINATION FIGURE;
    if d["recombination"]:
        a = []
        b = []
        for x in range(-50, 50, 1):
            y = x ** 2 + 2 * x + 2
            a.append(x)
            b.append(y)

        ax_symbol = fig.add_subplot(332)
        b = np.array(b)
        d = 500
        ax_symbol.plot(b - d, a, color='gray')
        ax_symbol.plot(-b + d, a, color='brown')
        ax_symbol.axis("off")

    plt.title("")

    return fig


def Execute(options):
    PrimerFilePath = os.path.join(options.inputDirectory, "PrimerData.csv")
    MatchedPrimerFilePath = os.path.join(options.inputDirectory, "MatchedPrimers.csv")
    PWMFilePath = os.path.join(options.inputDirectory, "PWMAnalysis.csv")

    # LOAD RELEVANT DATABASES;
    global PrimerData
    PrimerData = pd.read_csv(PrimerFilePath)
    PWMData = pd.read_csv(PWMFilePath)

    global MatchData
    MatchData = pd.read_csv(MatchedPrimerFilePath)

    # FETCH ORIGINAL HEATMAP GENOME LABELS;
    heatmapLabelsFilePath = os.path.join(
        options.inputDirectory,
        "heatmap_labels.npy"
    )

    global heatmapLabels
    heatmapLabels = np.load(heatmapLabelsFilePath)

    last = None

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


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-d",
                      dest="inputDirectory")

    options, args = parser.parse_args()

    Execute(options)
