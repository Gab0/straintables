#!/bin/python
import re
import os
import pandas as pd
import numpy as np
import array

import json
import subprocess

from optparse import OptionParser

import matplotlib.pyplot as plt

import fastcluster
import scipy

import detectMutations
import graphic

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GdkPixbuf

from matplotlib.figure import Figure
from matplotlib.backends.backend_gtk3agg import FigureCanvas
from matplotlib.backends.backend_gtk3 import (
    NavigationToolbar2GTK3 as NavigationToolbar)


class alignmentData():
    def __init__(self, inputDirectory):

        self.dataKeys = ["Unnamed: 0", "Unnamed: 1"]
        self.loadDataFiles(inputDirectory)
        self.inputDirectory = inputDirectory


    def loadDataFiles(self, inputDirectory):
        PrimerFilePath = os.path.join(inputDirectory, "PrimerData.csv")
        MatchedPrimerFilePath = os.path.join(inputDirectory, "MatchedPrimers.csv")
        PWMFilePath = os.path.join(inputDirectory, "PWMAnalysis.csv")

        # LOAD RELEVANT DATABASES;
        self.PrimerData = pd.read_csv(PrimerFilePath)
        self.PWMData = pd.read_csv(PWMFilePath)
        self.MatchData = pd.read_csv(MatchedPrimerFilePath)

        # FETCH ORIGINAL HEATMAP GENOME LABELS;
        heatmapLabelsFilePath = os.path.join(
            inputDirectory,
            "heatmap_labels.npy"
        )

        self.heatmapLabels = np.load(heatmapLabelsFilePath)
        # FETCH VIEWABLE DATA INDEXES;
        last = None
        self.allowedIndexes = []
        for I in range(self.PWMData.shape[0]):
            d = self.PWMData.iloc[I]
            a = d[self.dataKeys[0]]
            if a == last:
                continue
            self.allowedIndexes.append(I)
            last = a

    def findPWMDataRow(self, a, b):
        for k in range(self.PWMData.shape[0]):
            d = self.PWMData.iloc[k]
            names = [d[x] for x in self.dataKeys]
            if a in names and b in names:
                return d

        return None

    def buildArrayPath(self, f):
        return os.path.join(self.inputDirectory, f)

    def fetchLociList(self):
        for t in self.dataKeys:
            print(self.PWMData[t])
        allLoci = [list(self.PWMData[d]) for d in self.dataKeys]
        print(allLoci)
        allLoci = [j for s in allLoci for j in s]

        return list(set(allLoci))


# COMPLEX DISSIMILARITY MATRIX VIEWER GTK APPLICATION;
class matrixViewer():
    def __init__(self, inputDirectory=None):

        # INITIALIZE STATE VARIABLES;
        self.labelColorsOn = 1
        self.index = 0
        self.zoomedPlot = None
        self.swap = False

        self.infoText = Gtk.Label.new("")


        self.drawPlot = plotPwmIndex

        # INITIALIZE GTK WINDOW;
        win = Gtk.Window()
        win.connect("destroy", lambda x: Gtk.main_quit())
        win.set_default_size(400, 300)
        win.set_title("linkageMapper - Walk Chromosome Result")


        # PREPARE BUTTON ICON IMAGE;
        dna_icon_left = loadImage(graphic.dna_icon)
        dna_icon_right = loadImage(graphic.dna_icon)
        swap_icon = loadImage(graphic.swap)

        # INITIALIZE TOP MENU BAR;
        self.topMenubar = Gtk.MenuBar()

        # FILE dropdown;
        self.menuFile = Gtk.Menu()

        btn_menuFile = Gtk.MenuItem(label="File")
        btn_menuFile.set_submenu(self.menuFile)

        self.topMenubar.append(btn_menuFile)

        loadAlignment = Gtk.MenuItem(label="Load Alignment")
        loadAlignment.connect("activate", self.selectFolderPath)
        self.menuFile.append(loadAlignment)


        # NAVIGATION dropdown;
        self.menuNav = Gtk.Menu()
        btn_menuNav = Gtk.MenuItem(label="Navigation")
        btn_menuNav.set_submenu(self.menuNav)

        self.topMenubar.append(btn_menuNav)

        btn_jumpTo = Gtk.MenuItem(label="Jump to Loci")
        btn_jumpTo.connect("activate", self.selectLocusNames)
        self.menuNav.append(btn_jumpTo)



        # INITIALIZE LOCUS NAVIGATION BUTTONS;
        self.openSequenceLeft = Gtk.Button()
        self.openSequenceLeft.connect("clicked", lambda d: self.launchAlignViewer(0))
        self.openSequenceLeft.add(dna_icon_left)

        self.openSequenceRight = Gtk.Button()
        self.openSequenceRight.connect("clicked", lambda d: self.launchAlignViewer(1))
        self.openSequenceRight.add(dna_icon_right)

        self.btn_back = Gtk.Button(image=Gtk.Image(stock=Gtk.STOCK_GO_BACK))
        self.btn_back.connect("clicked", self.nav_back)

        self.btn_next = Gtk.Button(image=Gtk.Image(stock=Gtk.STOCK_GO_FORWARD))
        self.btn_next.connect("clicked", self.nav_forward)

        self.btn_invert = Gtk.Button()
        self.btn_invert.connect("clicked", self.swapPlot)
        self.btn_invert.add(swap_icon)

        toggleColor = Gtk.ToggleButton(label=None, image=Gtk.Image(stock=Gtk.STOCK_COLOR_PICKER))
        toggleColor.set_tooltip_text("Show Matrix Label Colors")
        toggleColor.set_active(True)
        toggleColor.connect("clicked", self.toggleColor)

        # INITIALIZE PLOT FIGURE;
        self.figure = plt.figure()

        #self.drawPlot(self.figure, self.PWMData, 0)
        self.figurecanvas = FigureCanvas(self.figure)

        self.figurecanvas.mpl_connect('button_press_event', self.onclickCanvas)

        # BUILD INTERFACE;
        vbox = Gtk.VBox()
        vbox.pack_start(self.topMenubar, expand=False, fill=False, padding=0)

        vbox.pack_start(self.figurecanvas, expand=True, fill=True, padding=0)

        # SHOW LOCUS NAVIGATION TOOLBAR;
        buttonBox = Gtk.Grid()

        self.btn_back.set_hexpand(True)
        self.btn_next.set_hexpand(True)

        buttonBox.add(self.openSequenceLeft)
        buttonBox.attach(self.btn_back, 1, 0, 3, 1)
        buttonBox.attach(self.btn_invert, 4, 0, 2, 1)
        buttonBox.attach(self.btn_next, 6, 0, 3, 1)
        buttonBox.attach(self.openSequenceRight, 9, 0, 2, 1)

        vbox.pack_start(buttonBox, expand=False, fill=True, padding=0)

        # MODIFY MATPLOTLIB TOOLBAR;
        self.toolbar = NavigationToolbar(self.figurecanvas, win)
        self.toolbar.set_history_buttons()

        # SET BOTTOM TOOLBAR, WHICH INCLUDE MATPLOTLIB BAR;
        panelBox = Gtk.HBox(homogeneous=False, spacing=2)


        panelBox.pack_start(self.toolbar, expand=False, fill=False, padding=0)
        panelBox.pack_start(toggleColor, expand=False, fill=False, padding=0)

        panelBox.add(self.infoText)

        vbox.pack_start(panelBox, False, False, 0)

        # SHOW ALL;
        win.add(vbox)
        win.show_all()

        # HIDE THIS ANNOYING THING.
        self.toolbar.message.hide()

        self.loadNewFolder(inputDirectory)
        # LAUNCH!
        Gtk.main()

    def loadNewFolder(self, inputDirectory):
        if inputDirectory:
            self.alnData = alignmentData(inputDirectory)
            self.infoText.set_text(self.alnData.inputDirectory)
            self.changeView(self.index)
        else:
            self.alnData = None

    def cycleIndexes(self, amt):
        if not self.alnData:
            return

        self.index += amt

        while self.index not in self.alnData.allowedIndexes:
            self.index += amt
            if self.index < 0:
                self.index = self.alnData.PWMData.shape[0] - 1
            if self.index > max(self.alnData.allowedIndexes):
                self.index = min(self.alnData.allowedIndexes)

    def nav_forward(self, d):
        self.swap = False
        self.cycleIndexes(1)
        self.changeView(self.index)

    def nav_back(self, d):
        self.swap = False
        self.cycleIndexes(-1)
        self.changeView(self.index)

    def toggleColor(self, d):
        self.labelColorsOn = 1 - self.labelColorsOn
        self.changeView(self.index)

    def swapPlot(self, d):
        self.swap = 1 - self.swap
        self.changeView(self.index)

    def changeView(self, I):
        self.figure.clf()
        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.drawPlot(self.figure, self.alnData, a, b, swap=self.swap, showLabelColors=self.labelColorsOn)
            try:
                self.figurecanvas.draw()
                self.figurecanvas.flush_events()
            except AttributeError:
                pass

            LocusNames = self.getLocusNames()

            self.openSequenceLeft.set_tooltip_text("View Alignment For %s" % LocusNames[0])
            self.openSequenceRight.set_tooltip_text("View Alignment For %s" % LocusNames[1])

    def onclickCanvas(self, event):
        # print(event)
        if event.inaxes:
            Axis = event.inaxes
            #if event.button is 1:
            if self.zoomedPlot is None:
                Axis._orig_position = Axis.get_position()
                Axis.set_position([0, 0, 1, 1])
                self.zoomedPlot = Axis

                if False:
                    self.figure.clf()
                    print(Axis)
                    self.figure.axes.append(Axis)
                    print(self.figure.axes)
                    self.figurecanvas.draw()
                    #self.figurecanvas.flush_events()


                for otherAxis in event.canvas.figure.axes:
                    if otherAxis is not Axis:
                        otherAxis.set_visible(False)

            #if event.button is 3:

            elif self.zoomedPlot is not None:
                # JUST REDRAW... SLOWER BUT GUARANTEED (matplotlib is mystical);
                self.changeView(self.index)
                self.zoomedPlot = None
                """
                    self.zoomedPlot.set_position(self.zoomedPlot._orig_position)
                    for Axis in event.canvas.figure.axes:
                        Axis.set_visible(True)
                    self.zoomedPlot = None
                    """
        else:
            print("OFF AXIS.;")
        """
            bbox = self.figure.get_window_extent().transformed(
                self.figure.dpi_scale_trans.inverted())
            width, height = bbox.width*self.figure.dpi, bbox.height*self.figure.dpi

            side = 0 if event.x < width // 2 else 1

            #self.launchAlignViewer(side)
        """

    def getLocusNames(self, fullName=False):
        Data = self.alnData.PWMData.iloc[self.index]
        KeyNames = ["Unnamed: 0", "Unnamed: 1"]
        locusNames = [Data[kn] for kn in KeyNames]
        if not fullName:
            locusNames = [n.replace(".npy", "") for n in locusNames]
        return locusNames


    def launchAlignViewer(self, side):
        LocusNames = self.getLocusNames()
        LocusName = LocusNames[side]

        alignmentFilePath = os.path.join(options.inputDirectory, LocusName)
        command = ["aliview", alignmentFilePath]
        print(command)

        subprocess.run(command)

    def selectFolderPath(self, n):
        def onResponse(widget, response):

            if response:
                inputDirectory = widget.get_filename()
                self.loadNewFolder(inputDirectory)
            widget.destroy()

        a = Gtk.FileChooserDialog(
            title="Select Results Folder",
            parent=None,
            action=Gtk.FileChooserAction.SELECT_FOLDER)

        a.add_button(Gtk.STOCK_OPEN, response_id=1)
        a.add_button(Gtk.STOCK_CANCEL, response_id=0)

        a.connect("response", onResponse)
        print(a)
        a.show()

    def selectLocusNames(self, n):
        if not self.alnData:
            return

        def onResponse(widget, response):

            print(dir(widget.left_choice))
            if response:
                a = widget.allLoci[widget.left_choice.get_active()]
                b = widget.allLoci[widget.right_choice.get_active()]
                self.figure.clf()
                self.drawPlot(self.figure, self.alnData, a, b)

            widget.destroy()

        a = Gtk.Dialog(
            title="Select Locus Names",
            parent=None)

        layout = Gtk.Grid()
        a.allLoci = self.alnData.fetchLociList()
        print(a.allLoci)

        optsAllLoci = Gtk.ListStore(str)

        for l in a.allLoci:
            optsAllLoci.append([l])


        a.left_choice = Gtk.ComboBox.new_with_model_and_entry(optsAllLoci)

        a.right_choice = Gtk.ComboBox.new_with_model_and_entry(optsAllLoci)

        a.left_choice.set_entry_text_column(0)
        a.right_choice.set_entry_text_column(0)
        layout.attach(a.left_choice, 0, 0, 1, 1)
        layout.attach(a.right_choice,0, 1, 1, 1)
        a.add(layout)
        box = a.get_content_area()
        box.add(layout)
        a.add_button(Gtk.STOCK_APPLY, response_id=1)
        a.add_button(Gtk.STOCK_CANCEL, response_id=0)
        a.set_default_size(150, 100)

        a.connect("response", onResponse)
        a.show_all()

def loadImage(source_image):
    bsource_image = source_image.tobytes()
    dnd = array.array('B', bsource_image)
    width, height = source_image.size
    image_pixelbuffer = GdkPixbuf.Pixbuf.new_from_data(dnd, GdkPixbuf.Colorspace.RGB, True, 8, width, height, width * 4)

    output_image = Gtk.Image()
    output_image.set_from_pixbuf(image_pixelbuffer)
    return output_image




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
        left = int(Z[cur_index - N, 0])
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


def singleLocusStatus(alnData, axis, locus_name):

    # FETCH HEALTH SCORE FOR LOCUS;
    locus_identifier = locus_name.replace("LOCI_", "")
    Health = alnData.MatchData[alnData.MatchData.LocusName == locus_identifier]
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


# well... this function can be..... simplified.
def matchPairOfClusterOutputData(clusterOutputData, Verbose=True):
    masterData = clusterOutputData[0]
    slaveData = clusterOutputData[1]

    if Verbose:
        print("Reordering Cluster Groups...")
        print()
        print("Original:")
        print(json.dumps(clusterOutputData, indent=2))
        print()

    # STEP ZERO: REORDER BY GROUP SIZE;
    def byGroupSize(Data):
        _Data = [(k, v) for k, v in Data.items()]
        _Data = sorted(_Data, key=lambda x: len(x[1]), reverse=True)
        _Data = [(idx, value[1]) for idx, value in enumerate(_Data)]
        _Data = dict(_Data)

        return _Data

    masterData = byGroupSize(masterData)
    slaveData = byGroupSize(slaveData)

    # 1st STEP: COMPUTE SCORES;
    def computeScores(masterData, slaveData):
        keyScores = [[0 for y in slaveData if len(slaveData[y]) > 1]
                     for x in masterData if len(masterData[x]) > 1]
        for mkey in range(len(keyScores)):
            for skey in range(len(keyScores[mkey])):
                score = 0
                for ind in masterData[mkey]:
                    if ind in slaveData[skey]:
                        score += 1
                #score = score / len(masterData[mkey])
                keyScores[mkey][skey] = score

        return keyScores

    keyScores = computeScores(masterData, slaveData)

    # DEBUG: VIEW SCORES;
    # clusterOutputData = [masterData, slaveData]
    # print(json.dumps(clusterOutputData, indent=2))
    for i in keyScores:
        for j in i:
            print("%i " % j, end="")
        print()

    print(keyScores)

    # 2nd STEP: SWAP POSITIONS @ SLAVE;
    TargetReplacementCount = (1, 4)
    ReplacementCount = 0
    Replaced = []
    for k in range(TargetReplacementCount[1]):
        for mkey in range(len(keyScores)):
            if mkey in Replaced:
                continue

            MAX = max(keyScores[mkey])
            MIN = min(keyScores[mkey])
            if len(keyScores[mkey]) >= len(keyScores):
                if keyScores[mkey][mkey] == max(keyScores[mkey]):
                    print("GOOD!")
                    continue
            for skey in range(len(keyScores[mkey])):
                if keyScores[mkey][skey] == MAX:
                    try:
                        poorIndex = mkey
                        goodIndex = skey

                        good = slaveData[goodIndex]
                        poor = slaveData[poorIndex]

                        slaveData[poorIndex] = good
                        slaveData[goodIndex] = poor

                        keyScores = computeScores(masterData, slaveData)
                        print("Replaced! %i %i" % (mkey, skey))
                    except KeyError as e:
                        print(e)

                    ReplacementCount += 1
                    Replaced.append(mkey)
        ReplacementCount += 1
        if ReplacementCount >= TargetReplacementCount[0]:
            break

    # 3rd STEP: CREATE NEW GROUPS @ SLAVE TO NOT REPEAT GROUPS FROM MASTER;
    # IF THEY DON'T SHARE A SINGLE GENOME;

    def doShareCommon(masterDataEntry, slaveDataEntry):
        shareCommon = False
        for Individue in masterDataEntry:
            if Individue in slaveDataEntry:
                shareCommon = True
                break
        return shareCommon

    masterDataLength = len(list(masterData.keys()))
    slaveDataLength = len(list(slaveData.keys()))
    for k in range(min(masterDataLength, slaveDataLength)):
        if not len(masterData[k]) > 1 or k not in masterData.keys():
            continue
        if not len(slaveData[k]) > 1 or k not in slaveData.keys():
            continue

        shareCommon = doShareCommon(masterData[k], slaveData[k])

        # IF NO GENOME IS SHARED, THIS STEP WILL END WITH THE DELETION OF K KEY FROM SLAVE DATA;
        # AND SLAVE DATA [K] WILL BE REASSIGNED TO RIGHT BEFORE THE LONE INDIVIDUES ARE, OR TO THE LAST POSITION;
        if not shareCommon:
            # Reassign possible existing group key to new:
            reassignedKey = None
            lastResourceKey = max(list(slaveData.keys())) + 1
            for skey in range(lastResourceKey):
                if skey in slaveData.keys():
                    if len(slaveData[skey]) == 1:
                        reassignedKey = skey
                        slaveData[lastResourceKey] = slaveData[reassignedKey]
                        slaveData[reassignedKey] = slaveData[k]
                        break
                elif skey in masterData.keys():
                    if doShareCommon(masterData[skey], slaveData[k]):
                        reassignedKey = skey
                        slaveData[skey] = slaveData[k]
                        break

            if reassignedKey is None:
                slaveData[lastResourceKey] = slaveData[k]

            if k in slaveData.keys():
                del slaveData[k]

    # clusterOutputData = [masterData, slaveData]
    # print(json.dumps(clusterOutputData, indent=2))

    #keyScores = computeScores(masterData, slaveData)
    # DEBUG: VIEW SCORES;
    clusterOutputData = [masterData, slaveData]
    print(json.dumps(clusterOutputData, indent=2))
    for i in keyScores:
        for j in i:
            print("%i " % j, end="")
        print()

    print(keyScores)

    return clusterOutputData


def plotPwmIndex(fig, alnData, a, b, swap=False, showLabelColors=True):

    if swap:
        c = b
        b = a
        a = c

    currentPWMData = alnData.findPWMDataRow(a, b)

    # walk loci by loci mode.

    # EXTRACR LOCUS NAMES;
    a_name, b_name = fixArrayFilename(a), fixArrayFilename(b)

    try:
        data = [
            alnData.PrimerData[alnData.PrimerData.Locus == name.replace("LOCI_", "")].iloc[0]
            for name in [a_name, b_name]
        ]
    except IndexError:
        print("Failure on %s" % a_name)

    # LOAD MATRIX DATA;
    ma = np.load(alnData.buildArrayPath(a))
    mb = np.load(alnData.buildArrayPath(b))

    # REORDERED MATRIXES;
    ordered_ma, matrix_order, B = compute_serial_matrix(ma, method="complete")
    ordered_mb = reorderMatrix(mb, matrix_order)
    orderedLabels = alnData.heatmapLabels[matrix_order]

    # plot;
    r_axis1 = createSubplot(fig, 331, a_name, ordered_ma, orderedLabels)
    r_axis2 = createSubplot(fig, 333, b_name, ordered_mb, orderedLabels)

    reordered_axis = [r_axis1, r_axis2]

    # ORIGINAL MATRIXES;
    # plot;
    o_axis1 = createSubplot(fig, 337, a_name, ma, alnData.heatmapLabels)
    o_axis2 = createSubplot(fig, 339, b_name, mb, alnData.heatmapLabels)

    original_axis = [o_axis1, o_axis2]

    # COLORIZE MATRIX LABELS BY MESHCLUSTER;
    if showLabelColors:
        # color map from matplotlib;
        colorMap = plt.get_cmap("tab20")

        GroupColors = [colorMap(x / 20)
                       for x in range(20)]


        # lower case greek letters for niceness;
        symbolMap = [chr(945 + x) for x in range(20)]

        clusterOutputData = [None for n in range(2)]
        # ITERATE LOCUS NAMES ON VIEW (TWO) iteration to load clusterOutputData;
        for N, LocusName in enumerate([a_name, b_name]):
            clusterFilePath = alnData.buildArrayPath(LocusName) + ".clst"
            if os.path.isfile(clusterFilePath):
                locusClusterOutputData = parseMeshcluster(clusterFilePath)
                clusterOutputData[N] = locusClusterOutputData

        # REORGANIZE CLUSTER OUTPUT DATA;
        if all(clusterOutputData):
            clusterOutputData = matchPairOfClusterOutputData(clusterOutputData)

        # NEW ITERATION OF LOCUS NAMES, TO APPLY CLUSTER OUTPUT DATA INTO VIEW;
        for N, LocusName in enumerate([a_name, b_name]):
            if clusterOutputData[N] is not None:
                NB_Groups = len(clusterOutputData[N].keys())
                for Axis in [reordered_axis[N], original_axis[N]]:

                    # COLORIZE LABELS;
                    axisLabels = list(zip(Axis.get_xticklabels(), Axis.get_yticklabels()))
                    for idx, (labelx, labely) in enumerate(axisLabels):
                        text = labelx.get_text()
                        for key in clusterOutputData[N].keys():
                            if text in clusterOutputData[N][key]:

                                # fetch current state of labels;
                                xcurrentState = Axis.get_xticklabels()
                                ycurrentState = Axis.get_yticklabels()

                                # BLACK COLOR AND NULL SYMBOL FOR ONE INDIVIDUAL GROUPS;
                                if len(clusterOutputData[N][key]) == 1:
                                    Symbol = " "
                                # COLOR AND GREEK SYMBOL FOR MULTI INDIVIDUAL GROUPS;
                                else:
                                    Symbol = symbolMap[key]
                                    labelx.set_color(GroupColors[key])
                                    labely.set_color(GroupColors[key])

                                # modify current label;
                                xcurrentState[idx] = Symbol + "   " + text
                                ycurrentState[idx] = text + "   " + Symbol

                                # apply new state to labels;
                                Axis.set_xticklabels(xcurrentState)
                                Axis.set_yticklabels(ycurrentState)
                                break

    # BUILD SHOWN INFO;
    distance = abs(data[0].PositionStart - data[1].PositionStart)

    INF_SYMBOL = chr(8734)
    Title = [
        "Distance = {:,} bp".format(distance),
        "%s vs %s" % (a_name, b_name),
        "Mantel=%.4f     p=%.4f" % (currentPWMData["mantel"], currentPWMData["mantel_p"]),
        "DIFF=%i" % currentPWMData["matrix_ranking_diff"],
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
    if "AlignmentHealth" in alnData.MatchData.keys():
        ax_ha = fig.add_subplot(334)
        ax_hb = fig.add_subplot(336)

        singleLocusStatus(alnData, ax_ha, a_name)
        singleLocusStatus(alnData, ax_hb, b_name)

        # Additional info on secondary axis DEPRECATED;
        if False:
            RecombinationMessage = "True" if currentPWMData["recombination"] else "False"
            Message = "Recombination? %s" % RecombinationMessage
            ax_hb.text(0.8, 1, s=Message)

    # RECOMBINATION FIGURE;
    if currentPWMData["recombination"]:
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

    # SHOW DATA;
    viewer = matrixViewer(options.inputDirectory)


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-d",
                      dest="inputDirectory")

    options, args = parser.parse_args()

    Execute(options)
