#!/bin/python

from Bio import SeqIO, AlignIO
from Bio import Data
from Bio.Align.Applications import ClustalwCommandline

from Bio import Align

import straintables.PrimerEngine.PrimerDesign as bfps
from straintables import OutputFile
from straintables.Database import genomeManager

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


def buildOutputName(GenomeName, RegionName, is_reverse, window):
    rev = "reverse_" if is_reverse else ""
    return "%s_%s_%s%i" % (GenomeName, RegionName, rev, window)


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


def runForWindow(options, protein, sequence, Window, Reverse):
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
    ID = buildOutputName(StrainName, region_name, Reverse, Window)
    DNA.id = ID
    DNA.description = ""

    PROT.id = ID
    PROT.description = ""

    # ProteinSequences.append(PROT)
    # DnaSequences.append(DNA)

    TestFilePrefix = "TEST_%s_%s" % (PROT.id, StrainName)

    OutDirectory = os.path.join(options.WorkingDirectory, "out_%s" % region_name)
    if not os.path.isdir(OutDirectory):
        os.mkdir(OutDirectory)

    SequencesAsStrings = [str(Sequence.seq) for Sequence in [PROT, protein]]

    alignscore = MakeTestAlignment(SequencesAsStrings)
    dndscore = MakeTestClustalAlignment(SequencesAsStrings,
                                        TestFilePrefix,
                                        region_name,
                                        OutDirectory)

    if region_name == sequence.id:
        print("check %s" % sequence.id)
        exit()

    if len(PROT.seq) > len(protein.seq):
        print("WARNING: protein fragment length is bigger than reference protein length!")

    print("%s: %s / %s" % (TestFilePrefix, dndscore, alignscore))

    return PROT, dndscore, alignscore


# NOT USED;
def ShowAlignment(TestFile):
    alan = subprocess.Popen(["alan", TestFile + ".aln"],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)

    res = alan.communicate()
    print(res[0])


def processAllTranslationWindows(options, protein, sequence):
    AlignmentScores = []
    for Reverse in range(2):
        for Window in range(3):
            WindowDescriptor = (Window, Reverse)
            (WindowSequence, dndscore, alignscore) = runForWindow(options,
                                                                  protein,
                                                                  sequence,
                                                                  Window,
                                                                  Reverse)
            AlignmentScores.append((WindowDescriptor, WindowSequence, dndscore))

    BestAlignment = sorted(AlignmentScores, key=lambda v: v[2])[0]
    return BestAlignment


def MakeTestAlignment(Sequences):

    # -- SETUP ALIGNER AND ITS SCORES;
    Aligner = Align.PairwiseAligner()

    Aligner.mode = "global"
    Aligner.open_gap_score = -1000
    Aligner.extend_gap_score = -1
    Aligner.match_score = 100
    # Aligner.gap_score = -100

    d = Aligner.align(*Sequences)

    return d.score


def MakeTestClustalAlignment(Sequences, TestFilePrefix, region_name, OutDirectory):

    TestFile = TestFilePrefix + ".fasta"
    TestFilePath = os.path.join(OutDirectory, TestFile)
    SeqIO.write(Sequences, open(TestFilePath, 'w'), format="fasta")

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


def AnalyzeRegion(options, RegionSequenceSource):

    region_name = options.RegionName

    if not region_name:
        print("Region name undefined.")
        exit(1)

    # p_name = "%s_prot.fasta" % region_name
    RegionSequencesFilename = "LOCI_%s.fasta" % region_name
    RegionSequencesFilepath = os.path.join(options.WorkingDirectory,
                                           RegionSequencesFilename)

    if not os.path.isfile(RegionSequencesFilepath):
        print("Region sequences file not found at %s." % RegionSequencesFilepath)

    source_seq = RegionSequenceSource.fetchGeneSequence(region_name)

    if source_seq is None:
        return 0

    protein = SeqIO.SeqRecord(source_seq.translate())

    protein.id = region_name
    protein.description = ""

    sequences = SeqIO.parse(RegionSequencesFilepath, format="fasta")
    SuccessSequences = 0
    TotalSequences = 0

    AllRegionSequences = []

    for sequence in sequences:
        TotalSequences += 1
        (RecommendedWindow, WindowSequence, score) =\
            processAllTranslationWindows(options, protein, sequence)
        AllRegionSequences.append(WindowSequence)

        print("Correct Window: %s" % RecommendedWindow)
        RecommendedWindow = None
        if score < 0.15:
            SuccessSequences += 1

    OutputProteinFilePath = os.path.join(options.WorkingDirectory,
                                         "Protein_%s.fasta" % region_name)

    with open(OutputProteinFilePath, 'w') as f:
        SeqIO.write(AllRegionSequences, f, format="fasta")
    cmd = ClustalwCommandline("clustalw2",
                              infile=OutputProteinFilePath
    )

    OutputProteinReferenceFilePath = os.path.join(options.WorkingDirectory,
                                                  "Protein_ref_%s.fasta" % region_name)

    with open(OutputProteinReferenceFilePath, 'w') as f:
        SeqIO.write(protein, f, format="fasta")

    cmd.seqnos = "ON"
    cmd()

    successPercentage = SuccessSequences / TotalSequences * 100
    print("Rate for %s: %.2f%%" % (region_name, successPercentage))
    return successPercentage


def runDirectory(options, RegionSequenceSource):
    WantedFileQuery = "LOCI_([\w\d]+).fasta"
    files = [
        f for f in os.listdir(options.WorkingDirectory)
        if re.findall(WantedFileQuery, f)
    ]

    for f in files:
        region_name = re.findall(WantedFileQuery, f)[0]

        opt = types.SimpleNamespace()

        opt.RegionName = region_name
        opt.WorkingDirectory = options.WorkingDirectory
        successPercentage = AnalyzeRegion(opt, RegionSequenceSource)
        if successPercentage < 100:
            with open("log", 'a') as f:
                f.write("%s with %.2f%%\n" % (region_name, successPercentage))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", dest="RegionName")
    parser.add_argument("-a", dest="RunDirectory", action="store_true")
    parser.add_argument("-d", dest="WorkingDirectory", default=".")
    options = parser.parse_args()
    return options


def Execute(options):
    InformationFile = OutputFile.AnalysisInformation(options.WorkingDirectory)
    InformationFile.read()

    AnnotationPath = InformationFile.content["annotation"]

    GenomeFilePaths = genomeManager.readGenomeFolder()

    GenomeFeatures = list(SeqIO.parse(AnnotationPath, format="genbank"))

    RegionSequenceSource = bfps.BruteForcePrimerSearcher(GenomeFeatures, GenomeFilePaths)

    if options.RunDirectory:
        runDirectory(options, RegionSequenceSource)
    else:
        AnalyzeRegion(options, RegionSequenceSource)


def main():
    options = parse_arguments()
    Execute(options)


if __name__ == "__main__":
    main()
