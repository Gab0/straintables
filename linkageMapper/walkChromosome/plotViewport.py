#!/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os

from . import matrixOperations, dissimilarityCluster

from linkageMapper import detectMutations


def fixArrayFilename(f):
    return f.split('.')[0]


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


def createMatrixSubplot(fig, position, name, matrix, labels):
    new_ax = fig.add_subplot(position)

    detectMutations.heatmapToAxis(matrix, new_ax, labels=labels)

    new_ax.set_xlabel(name)

    return new_ax


def cropLabels(Labels, maxSize=13, Replacer="..."):
    croppedLabels = []
    maxSize = maxSize - len(Replacer)
    for label in Labels:
        if len(label) > maxSize:
            crop_size = len(label) - maxSize
            crop_size += crop_size % 2
            crop_size //= 2

            mid_point = len(label) // 2

            allowed_side_size = mid_point - crop_size
            cropped = label[:allowed_side_size] + Replacer + label[-allowed_side_size:]
            
        else:
            cropped = label
        
        croppedLabels.append(cropped)
        
    return croppedLabels
    
        
def plotPwmIndex(fig, alnData, a, b, swap=False, showLabelColors=True):

    if swap:
        c = b
        b = a
        a = c

    currentPWMData = alnData.findPWMDataRow(a, b)

    # walk loci by loci mode.

    # EXTRACR LOCUS NAMES;
    a_name, b_name = fixArrayFilename(a), fixArrayFilename(b)
    print(a_name)

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
    ordered_ma, matrix_order, B = matrixOperations.compute_serial_matrix(ma, method="complete")
    ordered_mb = matrixOperations.reorderMatrix(mb, matrix_order)
    orderedLabels = alnData.heatmapLabels[matrix_order]
    
    # Crop label lengths;
    orderedLabels = cropLabels(orderedLabels)
    
    # plot;
    r_axis1 = createMatrixSubplot(fig, 331, a_name, ordered_ma, orderedLabels)
    r_axis2 = createMatrixSubplot(fig, 333, b_name, ordered_mb, orderedLabels)

    reordered_axis = [r_axis1, r_axis2]

    # ORIGINAL MATRIXES;
    # plot;
    o_axis1 = createMatrixSubplot(fig, 337, a_name, ma, alnData.heatmapLabels)
    o_axis2 = createMatrixSubplot(fig, 339, b_name, mb, alnData.heatmapLabels)

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
        abmatrix = [ma, mb]
        # ITERATE LOCUS NAMES ON VIEW (TWO) iteration to load clusterOutputData;
        for N, LocusName in enumerate([a_name, b_name]):
            clusterFilePath = alnData.buildArrayPath(LocusName) + ".clst"

            # MeShCluSt file exists.
            if os.path.isfile(clusterFilePath):
                locusClusterOutputData = dissimilarityCluster.parseMeshcluster(clusterFilePath)
            # Otherwise...
            else:
                locusClusterOutputData = dissimilarityCluster.fromDissimilarityMatrix(abmatrix[N], alnData.heatmapLabels)

            # Assign obtained clusters;
            clusterOutputData[N] = locusClusterOutputData

        # REORGANIZE CLUSTER OUTPUT DATA;
        if all(clusterOutputData):
            clusterOutputData = dissimilarityCluster.matchPairOfClusterOutputData(clusterOutputData)

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
    if currentPWMData is not None:
        color_green = (0.1, 0.8, 0.1)
        color_red = (0.8, 0.1, 0.1)

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
        # PWM[RECOMBINATION] IS DEPRECATED.
        # if currentPWMData["recombination"]:
        try:
            Recombination = dissimilarityCluster.checkRecombination(
                clusterOutputData,
                orderedLabels,
                Threshold=0.4)
        except Exception as e:
            print(clusterOutputData)
            Recombination = [False]
            print("WARNING: Recombination failure!")
            print(e)
            raise

        def plotRecombinationPanel(ax, baseIndex):
            x_values = np.linspace(0, 10, 100)

            pre = 0.7
            div = 2
            mul = 2.1
            plot_53 = [baseIndex + np.sin(pre + mul * x) / div for x in x_values]
            plot_35 = [baseIndex - np.sin(pre + mul * x) / div for x in x_values]

            ax.plot(x_values, plot_53, color=color_red)
            ax.plot(x_values, plot_35, color=color_green)

        if any(Recombination):
            a = []
            b = []
            for x in range(-50, 50, 1):
                y = x ** 2 + 2 * x + 2
                a.append(x)
                b.append(y)

            ax_recombination = fig.add_subplot(332)
            dm = list(range(len(orderedLabels)))

            # Reverse recombination array because matrix plot indexes and normal plot indexes are reversed.
            for r, rec in enumerate(reversed(Recombination)):
                if rec:
                    plotRecombinationPanel(ax_recombination, r)

            ax_recombination.scatter([0 for x in dm], dm, color=color_green)
            ax_recombination.scatter([10 for x in dm], dm, color=color_red)

            b = np.array(b)
            d = 500
            #ax_symbol.plot(b - d, a, color='gray')
            #ax_symbol.plot(-b + d, a, color='brown')
            ax_recombination.axis("off")

    plt.title("")

    plt.subplots_adjust(top=0.79, bottom=0.03, left=0.06, right=1.00)

    return fig

