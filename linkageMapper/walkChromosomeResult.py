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

import walkChromosome


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
        OnlySequence = False
        if OnlySequence:
            last = None
            self.allowedIndexes = []
            for I in range(self.PWMData.shape[0]):
                d = self.PWMData.iloc[I]
                a = d[self.dataKeys[1]]
                if a == last:
                    continue
                self.allowedIndexes.append(I)
                last = a
        else:
            self.allowedIndexes = list(range(self.PWMData.shape[0]))

        print("Allowed: %s" % self.allowedIndexes)

    def findPWMDataRow(self, a, b):
        def setLength(w):
            return len(list(set(w)))

        for k in range(self.PWMData.shape[0]):
            d = self.PWMData.iloc[k]
            names = [d[x] for x in self.dataKeys]

            if a in names and b in names and setLength(names) == setLength([a, b]):
                print(d)
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


class locusNamesSelectionMenu(Gtk.Grid):
    def __init__(self, matrixViewer):

        Gtk.Grid.__init__(self)
        self.matrixViewer = matrixViewer

        self.optsAllLoci = Gtk.ListStore(str)
        self.left_choice = Gtk.ComboBox.new_with_model_and_entry(self.optsAllLoci)
        self.right_choice = Gtk.ComboBox.new_with_model_and_entry(self.optsAllLoci)

        self.switchAutomaticDropdownLocusJump(Target=True)

        self.left_choice.set_entry_text_column(0)

        self.right_choice.set_entry_text_column(0)

        Refresh = Gtk.Button(label=".")
        Refresh.connect("clicked", self.Refresh)

        self.attach(Gtk.Label.new("<-"), 0, 0, 1, 1)
        self.attach(self.left_choice, 1, 0, 1, 1)
        self.attach(Refresh, 2, 0, 1, 1)
        self.attach(self.right_choice, 3, 0, 1, 1)
        self.attach(Gtk.Label.new("->"), 4, 0, 1, 1)

        self.allLoci = None
        self.show_all()

    def Refresh(self, n):
        print(self.matrixViewer.figure)

        a = self.allLoci[self.left_choice.get_active()]
        b = self.allLoci[self.right_choice.get_active()]
        self.matrixViewer.changeView(a, b)

        print("REFRESHED %s %s" % (a, b))

    def Update(self, n=None):
        if self.matrixViewer.alnData:
            self.allLoci = self.matrixViewer.alnData.fetchLociList()

            self.optsAllLoci.clear()
            for l in self.allLoci:
                self.optsAllLoci.append([l])

    def switchAutomaticDropdownLocusJump(self, Target=False):
        if Target:
            self.left_handler = self.left_choice.connect("changed", self.Refresh)
            self.right_handler = self.right_choice.connect("changed", self.Refresh)
        else:
            self.left_choice.disconnect(self.left_handler)
            self.right_choice.disconnect(self.right_handler)

    def updateInfo(self, a, b):
        if self.allLoci:
            # DISCONNECT COMBOBOX SIGNALS;
            self.switchAutomaticDropdownLocusJump(Target=False)

            self.left_choice.set_active(self.allLoci.index(a))
            self.right_choice.set_active(self.allLoci.index(b))

            # RECOONECT COMBOBOX SIGNALS;
            self.switchAutomaticDropdownLocusJump(Target=True)


class LocusMapBar(Gtk.Grid):
    def __init__(self, alnData=None):
        Gtk.Grid.__init__(self)

    def loadData(self, alnData):
        self.LocusNames = ['a', 'b']
        self.LocusButtons = []
        for n in self.LocusNames:
            button = Gtk.ss
        pass

    def updateView(self):
        pass

# COMPLEX DISSIMILARITY MATRIX VIEWER GTK APPLICATION;
class matrixViewer():
    def __init__(self, inputDirectory=None):

        # INITIALIZE STATE VARIABLES;
        self.labelColorsOn = 1
        self.index = 0
        self.zoomedPlot = None
        self.swap = False
        self.alnData = None
        self.infoText = Gtk.Label.new("")

        self.drawPlot = walkChromosome.plotViewport.plotPwmIndex

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

        self.locusNavigator = locusNamesSelectionMenu(self)

        # NAVIGATION dropdown;
        self.menuNav = Gtk.Menu()
        btn_menuNav = Gtk.MenuItem(label="Navigation")
        btn_menuNav.set_submenu(self.menuNav)

        self.topMenubar.append(btn_menuNav)

        btn_jumpTo = Gtk.MenuItem(label="Jump to Loci")
        btn_jumpTo.connect("activate", lambda e: locusNamesSelectionMenu(self))
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

        self.infoText.set_hexpand(True)
        panelBox.pack_start(self.infoText, expand=True, fill=True, padding=0)
        panelBox.pack_end(self.locusNavigator, False, False, 0)

        vbox.pack_start(panelBox, expand=False, fill=False, padding=0)

        # SHOW ALL;
        win.add(vbox)
        win.show_all()

        # HIDE THIS ANNOYING THING.
        self.toolbar.message.hide()

        self.loadNewFolder(inputDirectory)
        self.cycleIndexes(0)
        # LAUNCH!
        Gtk.main()

    def loadNewFolder(self, inputDirectory):
        if inputDirectory:
            self.alnData = alignmentData(inputDirectory)
            self.infoText.set_text(self.alnData.inputDirectory)
            self.nav_forward(None)
            self.locusNavigator.Update()
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

        print(self.index)

    def nav_forward(self, d):
        self.swap = False
        self.cycleIndexes(1)

        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.changeView(a, b)

    def nav_back(self, d):
        self.swap = False
        self.cycleIndexes(-1)
        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.changeView(a, b)

    def toggleColor(self, d):
        self.labelColorsOn = 1 - self.labelColorsOn
        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.changeView(a, b)

    def swapPlot(self, d):
        self.swap = 1 - self.swap
        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.changeView(a, b)

    def changeView(self, a, b):
        self.figure.clf()
        if self.alnData:
            self.locusNavigator.updateInfo(a, b)
            self.drawPlot(self.figure, self.alnData, a, b, swap=self.swap, showLabelColors=self.labelColorsOn)

            self.figurecanvas.draw()
            self.figurecanvas.flush_events()

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
                #self.changeView(self.)
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


def loadImage(source_image):
    bsource_image = source_image.tobytes()
    dnd = array.array('B', bsource_image)
    width, height = source_image.size
    image_pixelbuffer = GdkPixbuf.Pixbuf.new_from_data(dnd, GdkPixbuf.Colorspace.RGB, True, 8, width, height, width * 4)

    output_image = Gtk.Image()
    output_image.set_from_pixbuf(image_pixelbuffer)
    return output_image



def Execute(options):

    # SHOW DATA;
    viewer = matrixViewer(options.inputDirectory)


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-d",
                      dest="inputDirectory")

    options, args = parser.parse_args()

    Execute(options)
