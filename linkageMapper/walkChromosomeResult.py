#!/bin/python

import re
import os
import pandas as pd
import numpy as np
import array

import json
import subprocess

from argparse import ArgumentParser

import matplotlib.pyplot as plt

import fastcluster
import scipy

from linkageMapper import detectMutations
from . import graphic
from . import walkChromosome

import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, Gdk, GdkPixbuf
import cairo

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
        OnlySequence = True
        if OnlySequence:
            last = None
            self.allowedIndexes = []
            for I in range(self.PWMData.shape[0]):
                d = self.PWMData.iloc[I]
                a = d[self.dataKeys[0]]
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

        #w = Gtk.Button.new(Gdk.RGBA(100,100,100,100))
        #btn = ColoredButton()
        #btn.changeColor("red")
        btn = Gtk.Label.new("<-")
        self.attach(btn, 0, 0, 1, 1)
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


class LocusMapBar(Gtk.DrawingArea):
    def __init__(self):
        Gtk.DrawingArea.__init__(self)

        self.connect("draw", self.draw)
        self.set_size_request(200, 10)

        self.LocusNames = []
        self.Active = []
        self.show_all()

    def drawCircle(self, ctx, color=None):

        ctx.translate(self.circleSize, self.circleSize)

        ctx.new_path()
        ctx.arc(0, 0, self.circleSize, 0, 2 * 3.14)
        ctx.close_path()

        if color:
            ctx.set_source_rgba(*color, 1.0)
            ctx.fill()

        ctx.translate(-self.circleSize, -self.circleSize)

    def draw(self, da, ctx):
        # print("DRAWING %s" % self.Active)

        availableWidth = self.get_allocation().width

        Size = len(self.LocusNames)

        self.circleSize = max(min(availableWidth / (Size * 3), 12), 6)


        ctx.set_source_rgb(0, 0, 0)

        ctx.set_line_width(self.circleSize / 4)
        ctx.set_tolerance(0.1)

        # FIRST ROW;
        ctx.set_line_join(cairo.LINE_JOIN_ROUND)

        ctx.save()
        ctx.new_path()

        ctx.translate(self.circleSize, self.circleSize)

        nb_rows = 1
        for k, locus in enumerate(self.LocusNames):
            color = (0, 0, 0)

            if len(self.Active):
                if locus in self.Active[0]:
                    color = (0.1, 0.8, 0.1)
                elif locus in self.Active[1]:
                    color = (0.8, 0.1, 0.1)

            self.drawCircle(ctx, color)
            ctx.translate(3 * self.circleSize, 0)

            position_x = ctx.get_matrix()[4]

            if position_x > availableWidth:
                ctx.restore()
                ctx.save()
                ctx.translate(self.circleSize, self.circleSize + 3 * self.circleSize * nb_rows)
                nb_rows += 1

        print(nb_rows)
        self.set_size_request(200, 4 * self.circleSize * nb_rows)
        ctx.restore()

    def loadData(self, alnData):
        self.LocusNames = list(alnData.MatchData["LocusName"])
        print(self.LocusNames)


# COMPLEX DISSIMILARITY MATRIX VIEWER GTK APPLICATION;
class MatrixViewer():
    def __init__(self, inputDirectory=None):

        self.Online = False

        # INITIALIZE STATE VARIABLES;
        self.labelColorsOn = 1
        self.index = 0
        self.zoomedPlot = None
        self.swap = False
        self.alnData = None
        self.infoText = Gtk.Label.new("")

        # SELECT DRAWING FUNCTION;
        self.drawPlot = walkChromosome.plotViewport.plotPwmIndex

        # INITIALIZE GTK WINDOW;
        self.Window = Gtk.Window()
        self.Window.connect("destroy", lambda x: Gtk.main_quit())
        self.Window.set_default_size(400, 300)
        self.Window.set_title("linkageMapper - Walk Chromosome Result")

        # Children structures;
        self.locusMap = LocusMapBar()
        self.locusNavigator = locusNamesSelectionMenu(self)

        # LOAD DATA;
        self.loadNewFolder(inputDirectory)


        #self.top_menu()
        self.locus_navigator()
        Layout = self.build_interface()


        # -- THIS IS REQUIRED TO THE FIG.TIGHT_LAYOUT() TO SETTLE DOWN. A MYSTERY...

        self.Online = True
        # SHOW ALL;
        self.Window.add(Layout)
        self.Window.show_all()

        # HIDE THIS ANNOYING THING.
        self.toolbar.message.hide()

        self.navigate(0)
        self.navigate(0)
        if self.alnData:
            self.Window.connect("show", lambda x: self.navigate(0))

        # LAUNCH
        Gtk.main()

    def top_menu(self):
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
        btn_jumpTo.connect("activate", lambda e: locusNamesSelectionMenu(self))
        self.menuNav.append(btn_jumpTo)

    def locus_navigator(self):
        # PREPARE BUTTON ICON IMAGE;
        dna_icon_left = loadImage(graphic.dna_icon)
        dna_icon_right = loadImage(graphic.dna_icon)
        swap_icon = loadImage(graphic.swap)

        # INITIALIZE LOCUS NAVIGATION BUTTONS;
        self.openSequenceLeft = Gtk.Button()
        self.openSequenceLeft.connect("clicked", lambda d: self.launchAlignViewer(0))
        self.openSequenceLeft.add(dna_icon_left)

        self.openSequenceRight = Gtk.Button()
        self.openSequenceRight.connect("clicked", lambda d: self.launchAlignViewer(1))
        self.openSequenceRight.add(dna_icon_right)

        self.btn_back = Gtk.Button(image=Gtk.Image(stock=Gtk.STOCK_GO_BACK))
        self.btn_back.connect("clicked", lambda d: self.navigate(-1))

        self.btn_next = Gtk.Button(image=Gtk.Image(stock=Gtk.STOCK_GO_FORWARD))
        self.btn_next.connect("clicked", lambda d: self.navigate(1))

        self.btn_invert = Gtk.Button()
        self.btn_invert.connect("clicked", self.swapPlot)
        self.btn_invert.add(swap_icon)

        self.btn_toggleColor = Gtk.ToggleButton(label=None, image=Gtk.Image(stock=Gtk.STOCK_COLOR_PICKER))
        self.btn_toggleColor.set_tooltip_text("Show Matrix Label Colors")
        self.btn_toggleColor.set_active(True)
        self.btn_toggleColor.connect("clicked", self.toggleColor)


        # INITIALIZE PLOT FIGURE;
        self.figure = plt.figure()
        self.figurecanvas = FigureCanvas(self.figure)

        self.figurecanvas.mpl_connect('button_press_event', self.onclickCanvas)

    def build_interface(self):
        # BUILD INTERFACE;
        vbox = Gtk.VBox()
        # vbox.pack_start(self.topMenubar, expand=False, fill=False, padding=0)

        vbox.pack_start(self.locusMap, expand=False, fill=False, padding=0)

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
        self.toolbar = NavigationToolbar(self.figurecanvas, self.Window)
        self.toolbar.set_history_buttons()

        # SET BOTTOM TOOLBAR, WHICH INCLUDE MATPLOTLIB BAR;
        panelBox = Gtk.HBox(homogeneous=False, spacing=2)

        self.loadAlignment = Gtk.Button(image=Gtk.Image(stock=Gtk.STOCK_DND_MULTIPLE))
        self.loadAlignment.connect("clicked", self.selectFolderPath)
        panelBox.pack_start(self.loadAlignment, expand=False, fill=False, padding=3)

        panelBox.pack_start(self.toolbar, expand=False, fill=False, padding=0)
        panelBox.pack_start(self.btn_toggleColor, expand=False, fill=False, padding=0)

        self.infoText.set_hexpand(True)
        panelBox.pack_start(self.infoText, expand=True, fill=True, padding=0)
        panelBox.pack_end(self.locusNavigator, False, False, 0)

        vbox.pack_start(panelBox, expand=False, fill=False, padding=0)

        return vbox

    def loadNewFolder(self, inputDirectory):
        if inputDirectory:
            self.alnData = alignmentData(inputDirectory)
            self.infoText.set_text(self.alnData.inputDirectory)
            self.locusNavigator.Update()
            self.locusMap.loadData(self.alnData)
            if self.Online:
                self.navigate(0)
            self.inputDirectory = inputDirectory
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

    def navigate(self, value):
        self.swap = False
        self.cycleIndexes(value)
        if self.alnData:
            a, b = self.getLocusNames(fullName=True)
            self.changeView(a, b)

    def toggleColor(self, d):
        self.labelColorsOn = 1 - self.labelColorsOn
        self.navigate(0)

    def swapPlot(self, d):
        if self.alnData:
            self.swap = 1 - self.swap
            loci = self.getLocusNames(fullName=True)
            if self.swap:
                loci.reverse()
            self.changeView(*loci)

    def changeView(self, a, b):
        self.figure.clf()
        if self.alnData:
            # UPDATE LOCUS NAVIGATOR;
            self.locusNavigator.updateInfo(a, b)

            # UPDATE LOCUS MAP;
            self.locusMap.Active = [a, b]
            self.locusMap.queue_draw()

            self.drawPlot(self.figure, self.alnData, a, b, showLabelColors=self.labelColorsOn)

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

        alignmentFilePath = os.path.join(self.inputDirectory, LocusName)
        command = ["aliview", alignmentFilePath]
        print(command)

        subprocess.run(command)

    def selectFolderPath(self, n):
        a = folderPathSelector(self)


class folderPathSelector(Gtk.FileChooserDialog):
    def onResponse(self, widget, response):

        if response:
            inputDirectory = widget.get_filename()
            try:
                self.mainApplication.loadNewFolder(inputDirectory)
                widget.destroy()
            except Exception as e:
                ERR = e
                widget.errorMessage.set_text("Error: %s" % e)
        else:
            widget.destroy()

    def __init__(self, mainApplication):
        Gtk.FileChooserDialog.__init__(self,
                                       title="Select Results Folder",
                                       parent=None,
                                       action=Gtk.FileChooserAction.SELECT_FOLDER)

        self.mainApplication = mainApplication
        self.errorMessage = Gtk.Label.new("Select file.")
        self.errorMessage.set_hexpand(True)
        ok = Gtk.Button(Gtk.STOCK_OPEN)
        cancel = Gtk.Button(Gtk.STOCK_CANCEL)

        buttonGrid = Gtk.Grid()
        buttonGrid.attach(self.errorMessage, 0, 0, 1, 1)
        #buttonGrid.attach(ok, 0, 1, 1, 1)
        #buttonGrid.attach(cancel, 0, 2, 1, 1)

        self.vbox.pack_start(buttonGrid, False, False, 0)

        self.add_button(Gtk.STOCK_OPEN, response_id=1)
        self.add_button(Gtk.STOCK_CANCEL, response_id=0)
        self.connect("response", self.onResponse)

        self.show_all()


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
    if os.path.isfile(os.path.join(options.inputDirectory, "PrimerData.csv")):
        viewer = MatrixViewer(options.inputDirectory)
    else:
        print("Analysis file not found at input directory %s." % options.inputDirectory)


def parse_args():
    parser = ArgumentParser()

    parser.add_argument('inputDir', type=str, nargs=1, metavar="inputDirectory",
                        help='inputDirectory')

    parser.add_argument("-d", metavar="inputDirectory",
                        dest="inputDirectory")

    options = parser.parse_args()

    if not options.inputDirectory:
        options.inputDirectory = options.inputDir[0]

    return options


def main():
    options = parse_args()
    Execute(options)


if __name__ == "__main__":
    main()
