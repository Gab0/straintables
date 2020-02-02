#!/bin/python

import numpy as np


class LabelGroup():
    def __init__(self, baseNames):
        self.base = baseNames
        self.cropped = self.crop(self.base)

        # lowercase greek letters for niceness;
        self.clusterSymbolMap = [chr(945 + x) for x in range(55)]
        print(self.cropped)

    @staticmethod
    def crop(Labels, maxSize=13, Replacer="..."):
        croppedLabels = []
        maxSize -= len(Replacer)
        for label in Labels:
            if len(label) > maxSize:
                crop_size = len(label) - maxSize
                crop_size += crop_size % 2
                crop_size //= 2

                mid_point = len(label) // 2

                allowed_side_size = mid_point - crop_size

                cropped = label[:allowed_side_size]
                cropped += Replacer
                cropped += label[-allowed_side_size:]

            else:
                cropped = label

            croppedLabels.append(cropped)

        return croppedLabels

    def clusterize(self, clusterGroup):
        Cluster = [None for z in self.base]
        for n in clusterGroup.keys():
            if len(clusterGroup[n]) > 1:
                for member in clusterGroup[n]:
                    idx = None
                    for l, label in enumerate(self.base):
                        if label == member or label == member[:30]:
                            idx = l
                    if idx is not None:
                        Cluster[idx] = n

        return Cluster

    def get_labels(self, Cluster=[], symbolSide=0):
        Output = []

        symbolSideFormat = [
            "{symbol}{spacer}{label}",
            "{label}{spacer}{symbol}"
        ]

        for k, label in enumerate(self.cropped):
            if Cluster and Cluster[k] is not None:
                symbol = self.clusterSymbolMap[Cluster[k]]
            else:
                symbol = " "

            label_content = {
                'label': label,
                'spacer': " " * (15 - len(label)),
                'symbol': symbol
            }

            output_label = symbolSideFormat[symbolSide].format(**label_content)
            Output.append(output_label)

        return Output

    def get_ordered(self, reorderIndexes, **kwargs):
        r = np.array(self.get_labels(**kwargs))
        r = r[reorderIndexes]
        return list(r)
