#!/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os
import math

from . import matrixOperations, dissimilarityCluster
from . import MatrixPlot, MatrixLabelGroup


class PlotableMatrix():
    def __init__(self,
                 name,
                 matrix,
                 xlabels,
                 ylabels,
                 cluster,
                 colorizeSubplot,
                 MatrixParameters={}):
        self.name = name
        self.matrix = matrix
        self.xlabels = xlabels
        self.ylabels = ylabels
        self.cluster = cluster
        self.colorizeSubplot = colorizeSubplot
        self.MatrixParameters = MatrixParameters


def fixArrayFilename(f):
    return f.split('.')[0]


def reorderList(List, index_map):
    lst = np.array(List)[index_map]
    return list(lst)


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


def createMatrixSubplot(fig,
                        position,
                        plotableMatrix):
    new_ax = fig.add_subplot(position)

    MatrixPlot.drawMatrixOnAxis(
        fig,
        new_ax,
        plotableMatrix)

    return new_ax


def loadClusterPairData(alnData, a_name, b_name, abmatrix, Labels):

    clusterOutputData = [None for n in range(2)]
    # ITERATE LOCUS NAMES ON VIEW (TWO) iteration to load clusterOutputData;
    for N, region_name in enumerate([a_name, b_name]):
        clusterOutputData[N] = loadClusterData(alnData,
                                               region_name,
                                               abmatrix[N], Labels)

    # REORGANIZE CLUSTER OUTPUT DATA;
    if all(clusterOutputData):
        clusterOutputData =\
            dissimilarityCluster.matchPairOfClusterOutputData(
                clusterOutputData)

    return clusterOutputData


def loadClusterData(alnData, region_name, matrix, Labels):

    # Assign obtained clusters;
    clusterFilePath = alnData.buildArrayPath(region_name) + ".clst"

    # MeShCluSt file exists.
    if os.path.isfile(clusterFilePath):
        locusClusterOutputData =\
            dissimilarityCluster.parseMeshcluster(clusterFilePath)

    # Otherwise...
    else:
        locusClusterOutputData =\
            dissimilarityCluster.fromDissimilarityMatrix(matrix, Labels.base)

    return locusClusterOutputData


def makePlotCode(a, b, c):
    return a * 100 + b * 10 + c


def plotRegionBatch(fig,
                    alnData,
                    regionIndexes,
                    showLabelColors=True,
                    reorganizeIndex=None,
                    MatrixParameters={}):
    data = [
        alnData.MatchData["LocusName"].iloc[i]
        for i in regionIndexes
    ]

    Matrices = [np.load(alnData.buildArrayPath(a)) for a in data]

    Labels = MatrixLabelGroup.LabelGroup(alnData.heatmapLabels)
    Clusters = [
        loadClusterData(alnData, data[i], Matrices[i], Labels)
        for i in range(len(data))
    ]

    # Reorganize matrix logic;
    if reorganizeIndex is not None:
        guideData = alnData.MatchData["LocusName"].iloc[reorganizeIndex]
        guideMatrix = np.load(alnData.buildArrayPath(guideData))
        _, matrix_order, B =\
            matrixOperations.compute_serial_matrix(
                guideMatrix,
                method="complete"
            )

        Matrices = [
            matrixOperations.reorderMatrix(mat, matrix_order)
            for mat in Matrices
        ]

        Labels = MatrixLabelGroup.LabelGroup(
            alnData.heatmapLabels[matrix_order])

    AllAxis = []

    # Compute number of rows and columns for plot figure;
    NBL = len(data)
    NBROWS = min(2, math.ceil(NBL / 2))
    NBCOLS = math.ceil(NBL / NBROWS)

    AllPlots = []
    AllMatrices = []
    print("Plot Count: %i\nColumns: %i\nRows: %i" % (NBL, NBCOLS, NBROWS))
    for m, Matrix in enumerate(Matrices):
        print("Building...")
        PlotCode = makePlotCode(NBROWS, NBCOLS, m + 1)

        plotCluster = Labels.clusterize(Clusters[m])
        plotLabels = Labels.get_labels(Cluster=plotCluster)

        plot = fig.add_subplot(PlotCode)
        plotableMatrix = PlotableMatrix(data[m],
                                        Matrix,
                                        plotLabels,
                                        plotLabels,
                                        plotCluster,
                                        showLabelColors,
                                        MatrixParameters)
        MatrixPlot.drawMatrixOnAxis(fig,
                                    plot,
                                    plotableMatrix)

        AllMatrices.append(plotableMatrix)
        AllPlots.append(plot)
        AllAxis.append(plot)

    return AllMatrices


def MainDualRegionPlot(fig,
                       alnData,
                       regionIndexes,
                       showLabelColors=True,
                       MatrixParameters={}):

    # EXTRACR REGION NAMES;
    region_names = alnData.getRegionNamesFromIndex(regionIndexes)
    a_name, b_name = region_names

    currentPWMData = alnData.findPWMDataRow(*region_names)

    # LOAD MATRIX DATA;
    Matrices = [ma, mb] = [
        np.load(alnData.buildArrayPath(name))
        for name in region_names
    ]

    # Crop label lengths;
    Labels = MatrixLabelGroup.LabelGroup(alnData.heatmapLabels)

    ordered_ma, matrix_order, B =\
        matrixOperations.compute_serial_matrix(ma, method="complete")

    ordered_mb = matrixOperations.reorderMatrix(mb, matrix_order)

    OrderedMatrices = [ordered_ma, ordered_mb]

    # -- CLUSTER INFORMATION TO LABEL;
    abmatrix = [ma, mb]
    clusterOutputData = loadClusterPairData(
        alnData,
        a_name,
        b_name,
        abmatrix,
        Labels
    )

    LeftCluster = Labels.clusterize(clusterOutputData[0])
    RightCluster = Labels.clusterize(clusterOutputData[1])

    # -- DEFINE FONTSIZE FOR PLOT LABELS;
    if "fontsize" not in MatrixParameters.keys():
        MatrixParameters["fontsize"] = 40 / math.sqrt(ma.shape[0])

    AllMAT = []
    AllAxis = []
    for vertical in range(2):
        for horizontal in range(2):
            BaseCluster = [LeftCluster, RightCluster][horizontal]
            Name = region_names[horizontal]

            # -- Define matrix contents & cluster reorganization when required;
            if not vertical:
                CurrentMatrix = OrderedMatrices[horizontal]
                Cluster = reorderList(BaseCluster, matrix_order)
            else:
                CurrentMatrix = Matrices[horizontal]
                Cluster = BaseCluster

            # left plots have yticks on the right side.
            if horizontal == 0:
                MatrixParameters["YlabelsOnRight"] = True
            else:
                MatrixParameters["YlabelsOnRight"] = False

            Position = [vertical, horizontal]

            if Position == [0, 0]:
                xlabels = ylabels = Labels.get_ordered(
                    matrix_order,
                    Cluster=BaseCluster,
                    symbolSide=horizontal
                )

            elif Position == [0, 1]:
                xlabels = Labels.get_ordered(matrix_order,
                                             Cluster=BaseCluster, symbolSide=0)
                ylabels = Labels.get_ordered(matrix_order,
                                             Cluster=BaseCluster, symbolSide=1)

            elif Position == [1, 0]:
                xlabels = ylabels = Labels.get_labels(Cluster=BaseCluster)

            elif Position == [1, 1]:
                xlabels = Labels.get_labels(Cluster=BaseCluster, symbolSide=0)
                xlabels = Labels.get_labels(Cluster=BaseCluster, symbolSide=1)


            # -- Create Matrix and axes;
            MAT = PlotableMatrix(Name,
                                 CurrentMatrix,
                                 xlabels,
                                 ylabels,
                                 Cluster,
                                 showLabelColors,
                                 MatrixParameters)

            Code = 221 + 2 * vertical + horizontal
            Axis = createMatrixSubplot(fig, Code, MAT)


            AllMAT.append(MAT)
            AllAxis.append(Axis)

    # BUILD SHOWN INFO;
    if currentPWMData is not None:
        pass

    plt.title("")

    plt.subplots_adjust(top=0.79, bottom=0.03, left=0.06, right=1.00)
    # fig.tight_layout()

    return AllMAT


def plotRecombinationPanel(ax, baseIndex):

    color_green = (0.1, 0.8, 0.1)
    color_red = (0.8, 0.1, 0.1)
    x_values = np.linspace(0, 10, 100)

    pre = 0.7
    div = 2
    mul = 2.1

    plot_53 = [baseIndex + np.sin(pre + mul * x) / div
               for x in x_values]

    plot_35 = [baseIndex - np.sin(pre + mul * x) / div
               for x in x_values]

    ax.plot(x_values, plot_53, color=color_red)
    ax.plot(x_values, plot_35, color=color_green)


def RegionData(currentPWMData, MatchData, a_name, b_name):
    INF_SYMBOL = chr(8734)

    if MatchData[0]["Chromosome"] == MatchData[1]["Chromosome"]:
        try:
            distance = abs(MatchData[0]["PositionStart"] -
                           MatchData[1]["PositionStart"])
            distance = "{:,}".format(distance)
        except KeyError:
            print(MatchData)
            distance = INF_SYMBOL
    else:
        distance = INF_SYMBOL

    try:
        MantelData = "Mantel=%.4f     p=%.4f" % (currentPWMData["mantel"],
                                                 currentPWMData["mantel_p"])
        DiffData = "DIFF=%i" % currentPWMData["matrix_ranking_diff"]
    except TypeError:
        MantelData = "Mantel unknown."
        DiffData = "DIFF unknown."
    return [
        "Distance = %s bp" % distance,
        MantelData,
        DiffData,
        " "
    ]


def RegionInfoAxis(ax, Message):

    ax.text(
        0.2,
        0.6,
        s=Message,
        clip_on=False
    )

    ax.axis("off")


# DEPRECATED;
def AlignmentHealthAxis(ax_ha, ax_hb, alnData, currentPWMData, a_name, b_name):
    # ALIGNMENT HEALTH INFORMATION FIGURE;
    if "AlignmentHealth" in alnData.MatchData.keys():

        singleLocusStatus(alnData, ax_ha, a_name)
        singleLocusStatus(alnData, ax_hb, b_name)

        # Additional info on secondary axis DEPRECATED;
        """
            RecombinationMessage = "True" \
                if currentPWMData["recombination"] else "False"

            Message = "Recombination? %s" % RecombinationMessage
            ax_hb.text(0.8, 1, s=Message)
        """


# DEPRECATED;
def RecombinationAxis(fig, clusterOutputData, Labels, matrix_order):
    # RECOMBINATION FIGURE;

    color_green = (0.1, 0.8, 0.1)
    color_red = (0.8, 0.1, 0.1)
    try:
        Recombination = dissimilarityCluster.checkRecombination(
            clusterOutputData,
            Labels.get_ordered(matrix_order),
            Threshold=0.4)
    except Exception as e:
        print(clusterOutputData)
        # Recombination = [False]
        print("WARNING: Recombination failure!")
        print(e)
        raise

    if any(Recombination):
        a = []
        b = []
        for x in range(-50, 50, 1):
            y = x ** 2 + 2 * x + 2
            a.append(x)
            b.append(y)

        ax_recombination = fig.add_subplot(232)
        dm = list(range(len(Labels.base)))

        # Reverse recombination array,
        # because matrix plot indexes
        # and normal plot indexes are reversed.

        for r, rec in enumerate(reversed(Recombination)):
            if rec:
                plotRecombinationPanel(ax_recombination, r)

        ax_recombination.scatter([0 for x in dm], dm, color=color_green)
        ax_recombination.scatter([10 for x in dm], dm, color=color_red)

        ax_recombination.axis("off")
