#!/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os
import math

from . import matrixOperations, dissimilarityCluster
from . import MatrixPlot, MatrixLabelGroup


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
                        name,
                        matrix,
                        xlabels,
                        ylabels,
                        MatrixParameters={}):
    new_ax = fig.add_subplot(position)

    MatrixPlot.drawMatrixOnAxis(
        matrix,
        fig,
        new_ax,
        xlabels=xlabels,
        ylabels=ylabels,
        MatrixName=name,
        MatrixParameters=MatrixParameters
    )

    return new_ax


def colorizeSubplot(ax, Cluster):
    # color map from matplotlib;
    colorMap = plt.get_cmap("tab20")

    ClusterColors = [colorMap(x / 20)
                     for x in range(20)] * 4

    allLabels = enumerate(zip(ax.get_xticklabels(), ax.get_yticklabels()))
    for idx, (xlabel, ylabel) in allLabels:
        cluster = Cluster[idx]
        if cluster is not None:
            xlabel.set_color(ClusterColors[cluster])
            ylabel.set_color(ClusterColors[cluster])


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

    alignmentData = [
        alnData.AlignmentData[
            alnData.AlignmentData["LocusName"] == name].iloc[0]
        for name in data
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
    print("Plot Count: %i\nColumns: %i\nRows: %i" % (NBL, NBCOLS, NBROWS))
    for m, Matrix in enumerate(Matrices):
        print("Building...")
        PlotCode = makePlotCode(NBROWS, NBCOLS, m + 1)

        plotCluster = Labels.clusterize(Clusters[m])
        plotLabels = Labels.get_labels(Cluster=plotCluster)

        plot = fig.add_subplot(PlotCode)
        MatrixPlot.drawMatrixOnAxis(Matrix,
                                    fig,
                                    plot,
                                    xlabels=plotLabels,
                                    ylabels=plotLabels,
                                    MatrixName=data[m],
                                    MatrixParameters=MatrixParameters)

        AllPlots.append(plot)
        if showLabelColors:
            colorizeSubplot(plot, plotCluster)

        AllAxis.append(plot)

    for m, Matrix in enumerate(Matrices):
        axLabels = AllPlots[m].get_yticklabels()
        axLabel = axLabels[0]

        MatrixPlot.sequenceInfoOnAxis(
            AllPlots[m],
            reference=axLabel,
            nb_snp=alignmentData[m]["SNPCount"],
            aln_len=alignmentData[m]["AlignmentLength"],
            fontsize=MatrixParameters["fontsize"]
        )

    return fig


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
    [ma, mb] = [
        np.load(alnData.buildArrayPath(name))
        for name in region_names
    ]

    # Crop label lengths;
    Labels = MatrixLabelGroup.LabelGroup(alnData.heatmapLabels)

    ordered_ma, matrix_order, B =\
        matrixOperations.compute_serial_matrix(ma, method="complete")

    ordered_mb = matrixOperations.reorderMatrix(mb, matrix_order)

    # -- CLUSTER INFORMATION TO LABEL;
    abmatrix = [ma, mb]
    clusterOutputData = loadClusterPairData(alnData, a_name,
                                            b_name, abmatrix, Labels)

    LeftCluster = Labels.clusterize(clusterOutputData[0])
    RightCluster = Labels.clusterize(clusterOutputData[1])

    # -- DEFINE FONTSIZE FOR PLOT LABELS;
    if "fontsize" not in MatrixParameters.keys():
        MatrixParameters["fontsize"] = 40 / math.sqrt(ma.shape[0])

    # -- PLOT FOR REORDERED MATRICES (TOP);
    TA1_labels = Labels.get_ordered(matrix_order,
                                    Cluster=LeftCluster, symbolSide=0)
    top_axis1 = createMatrixSubplot(fig,
                                    221,
                                    a_name,
                                    ordered_ma,
                                    TA1_labels,
                                    TA1_labels,
                                    MatrixParameters=MatrixParameters)

    TA2_xlabels = Labels.get_ordered(matrix_order,
                                     Cluster=RightCluster, symbolSide=0)
    TA2_ylabels = Labels.get_ordered(matrix_order,
                                     Cluster=RightCluster, symbolSide=1)

    top_axis2 = createMatrixSubplot(fig,
                                    222,
                                    b_name,
                                    ordered_mb,
                                    TA2_xlabels,
                                    TA2_ylabels,
                                    MatrixParameters=MatrixParameters)

    # -- PLOT ORIGINAL MATRICES (BOTTOM);
    # plot;
    BA1_labels = Labels.get_labels(Cluster=LeftCluster)
    bottom_axis1 = createMatrixSubplot(fig,
                                       223,
                                       a_name,
                                       ma,
                                       BA1_labels,
                                       BA1_labels,
                                       MatrixParameters=MatrixParameters)

    BA2_xlabels = Labels.get_labels(Cluster=RightCluster, symbolSide=0)
    BA2_ylabels = Labels.get_labels(Cluster=RightCluster, symbolSide=1)
    bottom_axis2 = createMatrixSubplot(fig,
                                       224,
                                       b_name,
                                       mb,
                                       BA2_xlabels,
                                       BA2_ylabels,
                                       MatrixParameters=MatrixParameters)

    # left plots have yticks on the right side.
    top_axis1.yaxis.tick_right()
    bottom_axis1.yaxis.tick_right()

    # COLORIZE MATRIX LABELS BY MESHCLUSTER;
    if showLabelColors:
        colorizeSubplot(top_axis1, reorderList(LeftCluster, matrix_order))
        colorizeSubplot(bottom_axis1, LeftCluster)

        colorizeSubplot(top_axis2, reorderList(RightCluster, matrix_order))
        colorizeSubplot(bottom_axis2, RightCluster)

    # BUILD SHOWN INFO;
    if currentPWMData is not None:
        pass

    plt.title("")

    plt.subplots_adjust(top=0.79, bottom=0.03, left=0.06, right=1.00)
    # fig.tight_layout()

    return fig


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

    return [
        "Distance = %s bp" % distance,
        "Mantel=%.4f     p=%.4f" % (currentPWMData["mantel"],
                                    currentPWMData["mantel_p"]),
        "DIFF=%i" % currentPWMData["matrix_ranking_diff"],
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
