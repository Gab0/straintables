#!/bin/python

import os
import re
import math
import numpy as np
import matplotlib.cm
import matplotlib.pyplot as plt
from straintables import DrawGraphics, Definitions


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
        "fontsize": 9
    }
    # -- Process matrix parameters;
    for param in DefaultMatrixParameters.keys():
        if param not in MatrixParameters.keys():
            MatrixParameters[param] = DefaultMatrixParameters[param]

    return MatrixParameters


def drawMatrixOnAxis(MATRIX,
                     fig,
                     ax,
                     xlabels=None,
                     ylabels=None,
                     MatrixName=None,
                     MatrixParameters={}):

    if "fontsize" not in MatrixParameters.keys():
        MatrixParameters["fontsize"] = 40 / math.sqrt(MATRIX.shape[0])
    MatrixParameters = PreprocessMatrixParameters(MatrixParameters)

    if MatrixParameters["Normalize"]:
        MATRIX = normalizeMatrix(MATRIX, MatrixParameters)

    # -- Select colormap
    ColorMap = matplotlib.cm.get_cmap("binary")
    ax.matshow(MATRIX, cmap=ColorMap)

    SIZE = len(MATRIX)

    # MINOR TICKS -> GRID;
    DIV = SIZE // 3
    gridPoints = np.arange(0, SIZE, DIV)[1:-1] + 0.5

    ax.set_xticks(gridPoints, minor=True)
    ax.set_yticks(gridPoints, minor=True)

    # MAJOR TICKS -> LABELS;
    ax.set_xticks(range(SIZE))
    ax.set_yticks(range(SIZE))

    # SET LABELS;
    fontProperties = {
        'family': 'monospace',
        'fontsize': MatrixParameters["fontsize"]

    }

    # Calculate matrix cell size in pixels
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    cell_size = bbox.width
    fontProperties['fontsize'] = cell_size

    if xlabels is not None:
        ax.set_xticklabels(xlabels, fontProperties, rotation=90)
    if ylabels is not None:
        ax.set_yticklabels(ylabels, fontProperties)

    # Get approximate size for each cell.
    pos = ax.get_position()
    d = pos.height / SIZE
    print(d)

    if MatrixName:
        nameFontProperties = {}
        nameFontProperties.update(fontProperties)
        nameFontProperties['fontsize'] *= 4
        ax.set_xlabel(MatrixName, nameFontProperties)

    if MatrixParameters["showNumbers"]:
        valueFontProperties = fontProperties

        valueFontProperties['fontsize'] = cell_size / 2
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
