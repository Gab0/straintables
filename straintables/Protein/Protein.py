#!/bin/python

from Bio import SeqIO, AlignIO
from Bio import Data
from Bio.Align.Applications import ClustalwCommandline

from Bio import Align

import straintables.PrimerEngine.PrimerDesign as bfps

import copy
import argparse
import os
import subprocess
import re
import types


def parseDndFile(filepath, region_name):
    content = open(filepath).read()
    qr = "%s:([\d\.]+)\);" % region_name
    d = re.findall(qr, content)[0]
    return float(d)


def buildOutputName(Name, is_reverse, window):
    rev = "reverse_" if is_reverse else ""
    return "%s_%s%i" % (Name, rev, window)


class ReadFrameController():
    def __init__(self, Window, ReverseComplement):
        self.Window = Window
        self.ReverseComplement = ReverseComplement

    def apply(self, sequence):
        seq = copy.deepcopy(sequence)

        if self.ReverseComplement:
            seq = seq.reverse_complement()
        if self.Window > 0:
            seq = seq[self.Window:]

        TrimEnd = len(seq) % 3
        if TrimEnd:
            seq = seq[:-TrimEnd]
        assert(not len(seq) % 3)
        return seq


def runForWindow(protein, sequence, Window, Reverse):
    region_name = protein.id
    DNA = sequence

    RFC = ReadFrameController(Window, Reverse)
    DNA = RFC.apply(sequence)

    try:
        PROT = DNA.translate()
    except Data.CodonTable.TranslationError:
        print("TRANSLATION ERROR.")
        exit(1)

    StrainName = sequence.id
    DNA.id = buildOutputName(region_name, Reverse, Window)
    DNA.description = ""

    PROT.id = buildOutputName(region_name, Reverse, Window)
    PROT.description = ""

    # ProteinSequences.append(PROT)
    # DnaSequences.append(DNA)

    if False:
        # -- SETUP ALIGNER AND ITS SCORES;
        Aligner = Align.PairwiseAligner()

        Aligner.mode = "global"
        Aligner.open_gap_score = -1000
        Aligner.extend_gap_score = -1
        Aligner.match_score = 100
        # Aligner.gap_score = -100

        d = Aligner.align(PROT, protein)

        s = d.score

    OutDirectory = "out_%s" % region_name
    if not os.path.isdir(OutDirectory):
        os.mkdir(OutDirectory)

    if region_name == sequence.id:
        print("check %s" % sequence.id)
        exit()

    TestFilePrefix = "TEST_%s_%s" % (PROT.id, StrainName)
    TestFile = TestFilePrefix + ".fasta"

    TestFilePath = os.path.join(OutDirectory, TestFile)

    if True:
        ALIGN = [PROT, protein]
        if len(PROT) > len(protein):
            print("WARNING: protein fragment length is bigger than reference protein length!")

        SeqIO.write(ALIGN, open(TestFilePath, 'w'), format="fasta")

        Outfile = os.path.join(OutDirectory, TestFile + ".aln")
        cmd = ClustalwCommandline("clustalw2",
                                  infile=TestFilePath,
                                  outfile=Outfile)

        cmd.seqnos = "ON"
        cmd()

        dndfile = os.path.join(OutDirectory, TestFilePrefix + ".dnd")
        dndscore = parseDndFile(dndfile, region_name)
        if dndscore > 0.3:
            os.remove(TestFilePath)
        os.remove(dndfile)

        # x = AlignIO.read(Outfile, format='clustal')
        # print(str(x))

    print("%s: %s" % (TestFile, dndscore))

    if False:
        alan = subprocess.Popen(["alan", TestFile + ".aln"],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                shell=True)

        res = alan.communicate()
        print(res[0])

    return dndscore


def processAllTranslationWindows(protein, sequence):
    AlignmentScores = []
    for Reverse in range(2):
        for Window in range(3):
            WindowDescriptor = (Window, Reverse)
            dndscore = runForWindow(protein, sequence, Window, Reverse)
            AlignmentScores.append((WindowDescriptor, dndscore))

    BestAlignment = sorted(AlignmentScores, key=lambda v: v[1])[0]
    return BestAlignment


def AnalyzeRegion(options):

    ProteinSequences = []
    DnaSequences = []

    region_name = options.RegionName

    # p_name = "%s_prot.fasta" % region_name
    s_name = "LOCI_%s.fasta" % region_name

    if not os.path.isfile(s_name):
        print("Sequence file %s not found." % s_name)

    ResourceBasePath = "/home/gabs/linkage_analysis/toxoplasma"

    AnnotationPath = os.path.join(ResourceBasePath,
                                  "annotations/GCA_000006565.2_TGA4_genomic.gbff")
    GenomePaths = [os.path.join(ResourceBasePath, "genomes/ME49.fna")]

    GenomeFeatures = list(SeqIO.parse(AnnotationPath, format="genbank"))
    source = bfps.BruteForcePrimerSearcher(GenomeFeatures, GenomePaths)

    # protein = SeqIO.parse(p_name, format="fasta")
    # protein = list(protein)[0]

    source_seq = source.fetchGeneSequence(region_name, s_name)

    protein = SeqIO.SeqRecord(source_seq.translate())

    protein.id = region_name
    protein.description = ""

    sequences = SeqIO.parse(s_name, format="fasta")
    RecommendedWindow = None
    SuccessSequences = 0
    TotalSequences = 0
    for sequence in sequences:
        TotalSequences += 1
        if RecommendedWindow is None:
            RecommendedWindow, score = processAllTranslationWindows(protein, sequence)
            print(RecommendedWindow)
            RecommendedWindow = None
            if score < 0.15:
                SuccessSequences += 1

        else:
            runForWindow(sequence, *RecommendedWindow)

        if False:
            op_base = "OutputProteins_%s" % sequence.id
            os_base = "OutputSequences_%s" % sequence.id

            OutputProteins = op_base + ".fasta"
            OutputProteinsAlignment = op_base + ".aln"

            OutputSequences = os_base + ".fasta"
            OutputSequencesAlignment = os_base + ".aln"

            with open(OutputProteins, 'w') as f:
                SeqIO.write(ProteinSequences, f, format="fasta")
            with open(OutputSequences, 'w') as f:
                SeqIO.write(DnaSequences, f, format="fasta")

        if False:
            ClustalwCommandline("clustalw2",
                                infile=OutputProteins,
                                outfile=OutputProteinsAlignment)()
            ClustalwCommandline("clustalw2",
                                infile=OutputSequences,
                                outfile=OutputSequencesAlignment)()

    successPercentage = SuccessSequences / TotalSequences * 100
    print("Rate for %s: %.2f%%" % (region_name, successPercentage))
    return successPercentage


def runDirectory():
    files = [f for f in os.listdir() if f.endswith(".fasta")]

    for f in files:
        region_name = re.findall("_([\w\d]+).fasta", f)[0]
        opt = type('', (), {})()
        opt.RegionName = region_name
        successPercentage = AnalyzeRegion(opt)
        if successPercentage < 100:
            with open("log", 'a') as f:
                f.write("%s with %.2f%%\n" % (region_name, successPercentage))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", dest="RegionName")
    parser.add_argument("-a", dest="RunDirectory", action="store_true")
    options = parser.parse_args()
    return options


def main():
    options = parse_arguments()
    if options.RunDirectory:
        runDirectory()
    else:
        AnalyzeRegion(options)


if __name__ == "__main__":
    main()
