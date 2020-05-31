#!/bin/python

import os
import re
import math
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
from straintables import DrawGraphics, Definitions

matplotlib.use('ps')


# BUILD HEATMAP;
def createPdfHeatmap(MATRIX,
                     sequenceNames,
                     filename=None,
                     subtitle=None,
                     MatrixParameters={}):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    drawMatrixOnAxis(MATRIX,
                     ax,
                     xlabels=sequenceNames,
                     ylabels=sequenceNames,
                     MatrixParameters=MatrixParameters)

    # ax.grid(which='minor', color='r', linestyle='-', linewidth=2)

    if filename:
        plt.savefig(filename, bbox_inches='tight')

    FilenamePattern = r"%s(.+)\.aln" % Definitions.FastaRegionPrefix
    watermarkLabel = re.findall(FilenamePattern,
                                filename)
    # for loci similarity matrix;
    if watermarkLabel:
        watermarkLabel = watermarkLabel[0]
    # for other types of matrix;
    else:
        watermarkLabel = os.path.split(filename)[-1].split(".")[0]

    DrawGraphics.geneGraphs.watermarkAndSave(
        watermarkLabel,
        filename,
        subtitle=subtitle,
        verticalLabel=340
    )


def sequenceInfoOnAxis(ax, reference=None, nb_snp=0, aln_len=0, fontsize=10):
    Message = "# snp=%i\nlength=%i" % (nb_snp, aln_len)
    ax.annotate(
        Message,
        # xy=(-0.2, 1.2),
        xy=(0, 5),
        xycoords=reference,
        clip_on=False,
        ha='left',
        va='top',
        fontsize=fontsize
        )


def normalizeMatrix(MATRIX, parameters):

    MATRIX = MATRIX * parameters["pre_multiplier"]
    MODE = parameters["normalizer"]
    if MODE == 0:
        MATRIX = 1/(1+np.exp(-MATRIX))
    elif MODE == 1:
        MATRIX = np.tanh(MATRIX)
    elif MODE == 2:
        std = np.std(MATRIX)
        MATRIX = MATRIX * std
        MATRIX = np.tanh(MATRIX)

    return MATRIX


def PreprocessMatrixParameters(MatrixParameters):
    DefaultMatrixParameters = {
        "Normalize": False,
        "showNumbers": False,
        "fontsize": 9,
        "YlabelsOnRight": False,
    }
    # -- Process matrix parameters;
    for param in DefaultMatrixParameters.keys():
        if param not in MatrixParameters.keys():
            MatrixParameters[param] = DefaultMatrixParameters[param]

    return MatrixParameters


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


def drawMatrixOnAxis(fig,
                     ax,
                     plotableMatrix):

    MATRIX = plotableMatrix.matrix
    if "fontsize" not in plotableMatrix.MatrixParameters.keys():
        plotableMatrix.MatrixParameters["fontsize"] = 40 / math.sqrt(MATRIX.shape[0])
    MatrixParameters = PreprocessMatrixParameters(plotableMatrix.MatrixParameters)

    if MatrixParameters["Normalize"]:
        MATRIX = normalizeMatrix(MATRIX, MatrixParameters)

    # -- Select colormap
    ColorMap = matplotlib.cm.get_cmap("binary")
    ax.matshow(MATRIX, cmap=ColorMap)

    SIZE = len(MATRIX)

    # MINOR TICKS -> GRID;
    DIV = SIZE // 3
    gridPoints = np.arange(0, SIZE, DIV)[1:-1] + 0.5

    ax.set_xticks(gridPoints, minor=False)
    ax.set_yticks(gridPoints, minor=False)

    # MAJOR TICKS -> LABELS;
    ax.set_xticks(range(SIZE))
    ax.set_yticks(range(SIZE))

    if MatrixParameters["YlabelsOnRight"]:
        ax.yaxis.tick_right()

    # SET LABELS;
    fontProperties = {
        'family': 'monospace',
        'fontsize': MatrixParameters["fontsize"]

    }

    # Calculate matrix cell size in pixels
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    cell_size = bbox.width / SIZE ** 0.5 * 10
    print("Cell size: %.2f px" % cell_size)

    fontProperties['fontsize'] = min(fontProperties['fontsize'], cell_size)

    if plotableMatrix.xlabels is not None:
        ax.set_xticklabels(plotableMatrix.xlabels, fontProperties, rotation=90)
    if plotableMatrix.ylabels is not None:
        ax.set_yticklabels(plotableMatrix.ylabels, fontProperties)

    if plotableMatrix.name:
        nameFontProperties = {}
        nameFontProperties.update(fontProperties)
        nameFontProperties['fontsize'] *= 2
        ax.set_xlabel(plotableMatrix.name, nameFontProperties)

    if MatrixParameters["showNumbers"]:
        valueFontProperties = fontProperties

        valueFontProperties['fontsize'] = cell_size / 1.67
        # 2 * np.sqrt(fontProperties['fontsize'])
        Mean = np.mean(MATRIX)
        for i in range(MATRIX.shape[0]):
            for j in range(MATRIX.shape[1]):
                cell_value = MATRIX[i, j]
                txt_value = ("%.2f" % cell_value)[1:]

                # -- invert number colors for legibility
                if cell_value > Mean / 2:
                    color = 0
                else:
                    color = 255

                ax.text(j, i,
                        txt_value, valueFontProperties,
                        color=ColorMap(color), va='center', ha='center')

    ax.tick_params(
        axis='x',
        which='both',
        bottom=False,
        top=True,
        labelbottom=False
    )
    ax.tick_params(
        axis='both',
        which='minor',
        bottom=False,
        top=False,
        labelbottom=False
    )

    if plotableMatrix.colorizeSubplot:
        colorizeSubplot(ax, plotableMatrix.cluster)
