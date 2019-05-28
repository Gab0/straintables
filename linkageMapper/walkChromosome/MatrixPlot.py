#!/bin/python

import numpy as np


def heatmapToAxis(MATRIX, ax, xlabels=None, ylabels=None, fontsize=9, MatrixName=None):
    ax.matshow(MATRIX, cmap='binary')

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
        'fontsize': fontsize
    }
    if xlabels is not None:
        ax.set_xticklabels(xlabels, fontProperties, rotation=90)
    if ylabels is not None:
        ax.set_yticklabels(ylabels, fontProperties)

    if MatrixName:
        ax.set_xlabel(MatrixName, fontProperties)
