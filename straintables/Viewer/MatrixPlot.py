#!/bin/python

import numpy as np
import matplotlib.cm


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


def heatmapToAxis(MATRIX, ax, xlabels=None,
                  ylabels=None, fontsize=9,
                  MatrixName=None, MatrixParameters=None):

    # -- Process matrix parameters;
    if MatrixParameters is None:
        MatrixParameters = {
            "Normalize": False,
            "showNumbers": False
        }

    if MatrixParameters["Normalize"]:
        MATRIX = normalizeMatrix(MATRIX, MatrixParameters)

    ColorMap = matplotlib.cm.get_cmap("binary")
    ax.matshow(MATRIX, cmap=ColorMap)

    SIZE = len(MATRIX)
    print(MATRIX)
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
        'fontsize': fontsize
    }

    if xlabels is not None:
        ax.set_xticklabels(xlabels, fontProperties, rotation=90)
    if ylabels is not None:
        ax.set_yticklabels(ylabels, fontProperties)

    # Get approximate size for each cell.
    pos = ax.get_position()
    d = pos.height / MATRIX.shape[0]

    if MatrixName:
        ax.set_xlabel(MatrixName, fontProperties)

    if MatrixParameters["showNumbers"]:
        valueFontProperties = fontProperties

        valueFontProperties['fontsize'] = d * 430
        #2 * np.sqrt(fontProperties['fontsize'])
        Mean = np.mean(MATRIX)
        for i in range(MATRIX.shape[0]):
            for j in range(MATRIX.shape[1]):
                value = MATRIX[i, j]
                svalue = ("%.2f" % value)[1:]

                # -- invert number colors for legibility
                if value > Mean / 2:
                    color = 0
                else:
                    color = 255

                ax.text(j, i,
                        svalue, valueFontProperties,
                        color=ColorMap(color), va='center', ha='center')
