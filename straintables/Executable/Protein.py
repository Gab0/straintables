#!/bin/python

from Bio import Seq, SeqIO, SeqRecord
from Bio import Data
from Bio.Align.Applications import ClustalOmegaCommandline

import straintables.PrimerEngine.PrimerDesign as bfps
from straintables import OutputFile, Definitions
from straintables.Database import genomeManager

import copy
import argparse
import os
import subprocess
import re

import importlib
import random

Description = """

This command will translate multifasta sequences
inside a straintables' WorkingDirectory.

It will search thru the six translation windows to find the optimal.

"""

StandardAlignerParameters = {}


def buildOutputName(GenomeName, RegionName, RFC):
    return "%s_%s_%s" % (GenomeName, RegionName, RFC.FNSTR())


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

    def __str__(self):
        Orientation = ["5'", "3'"]
        if self.ReverseComplement:
            Orientation = reversed(Orientation)
        return "%s: +%i" % ("".join(Orientation), self.Window)

    def FNSTR(self):
        return str(self).replace(" ", "_").replace("'", "l")


def FindORF(sequence):
    return re.findall(r"M\w+\*", sequence)


def evaluateTranslationWindow(options,
                              TemplateProtein,
                              QuerySequence,
                              RFC,
                              Verbose,
                              AlignerParameters):

    region_name = TemplateProtein.id

    if Verbose:
        print("Evaluating window %s" % RFC)

    DNA = RFC.apply(QuerySequence)

    try:
        PROT = DNA.translate(table=1)
    except Data.CodonTable.TranslationError:
        print("TRANSLATION ERROR.")
        exit(1)

    StrainName = QuerySequence.id
    ID = buildOutputName(StrainName, region_name, RFC)

    DNA.id = ID
    DNA.description = ""

    PROT.id = ID
    PROT.description = ""

    # ProteinSequences.append(PROT)
    # DnaSequences.append(DNA)

    TestFilePrefix = "TEST_%s_%s" % (PROT.id, StrainName)

    OutDirectory = os.path.join(
        options.WorkingDirectory, "out_%s" % region_name)

    if options.WriteFiles:
        if not os.path.isdir(OutDirectory):
            os.mkdir(OutDirectory)

    # SequencesAsStrings = [str(s.seq) for s in [PROT, TemplateProtein]]

    # alignscore = MakeTestAlignment(SequencesAsStrings, AlignerParameters)
    alignscore = 0

    if options.OnlyOpenReadingFrame:
        AllORFs = FindORF(str(PROT.seq))

        ORFScores = []
        if TemplateProtein is not None:
            for ORF in AllORFs:
                if Verbose:
                    print("Searching ORF> %s" % ORF)
                AlignedSequences = ([str(TemplateProtein.seq), ORF])

                Alignment = MakeTestAlignment(AlignedSequences,
                                              AlignerParameters)

                if Verbose:
                    print()
                    print(Alignment[0])
                    print()

                    orf_alignscore = Alignment.score
                else:
                    orf_alignscore = MakeAlignment(AlignedSequences)
                if Verbose:
                    print("ORF Score = %.2f" % orf_alignscore)
                ORFScores.append(orf_alignscore)

            # -- SELECT THE BEST SCORING ORF FOR CURRENT WINDOW;
            ORF_SCORES = zip(AllORFs, ORFScores)
            BestORF, BestORFScore = sorted(
                ORF_SCORES,
                key=lambda os: os[1],
                reverse=True)[0]

            Result = SeqRecord.SeqRecord(Seq.Seq(BestORF), id=ID)
        else:
            raise(NotImplementedError)

    else:
        raise(NotImplementedError)

    if region_name == QuerySequence.id:
        print("check %s" % QuerySequence.id)
        exit()

    if Verbose:
        if len(PROT.seq) > len(TemplateProtein.seq):
            print("WARNING: protein fragment length > " +
                  "reference protein length!")

        print("%s:\n  aln: %s\norfaln: %s" %
              (TestFilePrefix, alignscore, orf_alignscore))
        print("\n\n")

    return Result, alignscore, BestORFScore


# NOT USED;
def ShowAlignment(TestFile):
    alan = subprocess.Popen(["alan", TestFile + ".aln"],
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE,
                            shell=True)

    res = alan.communicate()
    print(res[0])


def SortAlignment(ScoreData):
    WD, WS, align_score = ScoreData
    # EOT = len([x for x in WS if x == "*"])

    return align_score


def processAllTranslationWindows(options,
                                 TemplateProtein,
                                 QueryProtein,
                                 Verbose=0,
                                 AlignerParameters=StandardAlignerParameters):
    AlignmentScores = []
    for Reverse in range(2):
        for Window in range(3):
            WindowDescriptor = (Window, Reverse)

            RFC = ReadFrameController(Window, Reverse)
            (WindowSequence, dndscore, alignscore) = evaluateTranslationWindow(
                options,
                TemplateProtein,
                QueryProtein,
                RFC,
                Verbose=Verbose,
                AlignerParameters=AlignerParameters
            )

            if Verbose:
                print()
                print(WindowSequence.seq)
                print()

            AlignmentScores.append(
                (WindowDescriptor, WindowSequence, alignscore)
            )

    BestAlignment = sorted(AlignmentScores, key=SortAlignment, reverse=True)[0]
    return BestAlignment


def EvaluateSequenceByUnknownBase(options, sequence):
    UnknownBaseCount = len([base for base in sequence if base == "N"])
    if UnknownBaseCount >= options.DiscardImperfectSequenceThreshold:
        pass
    LongestUnknownBaseCluster = 0
    UKBC = 0
    for base in sequence:
        if base == "N":
            UKBC += 1
        else:
            LongestUnknownBaseCluster = max(UKBC,
                                            LongestUnknownBaseCluster)
            UKBC = 0

    if LongestUnknownBaseCluster >= options.DiscardImperfectSequenceThreshold:
        return False

    return True


def MakeAlignment(SequenceAsStrings: str):
    """

    This is used to find out which translation window sequence is most similar
    to a reference protein sequence.
    Bio.Align.PairwiseAligner somehow failed to achieve this.

    """
    MainSequence, QuerySequence = sorted(SequenceAsStrings,
                                         key=lambda s: -len(s))

    FragmentSize = 3
    FragmentCount = 50
    MatchedFragments = 0

    if len(QuerySequence) < FragmentSize:
        return 0

    for w in range(FragmentCount):
        PossibleIndexes = range(max(len(QuerySequence) - FragmentSize, 3))
        F = random.choice(PossibleIndexes)
        J = F + FragmentSize
        MatchedFragments += QuerySequence[F:J] in MainSequence

    return MatchedFragments


def MakeTestAlignment(SequencesAsStrings: str, AlignerParameters: dict):

    # -- SETUP ALIGNER AND ITS SCORES;
    Align = importlib.import_module("Bio.Align")
    Aligner = Align.PairwiseAligner()

    Aligner.mode = "local"

    Aligner.__dict__.update(AlignerParameters)

    d = Aligner.align(*SequencesAsStrings)

    return d


def LoadSequences(options, region_name):
    # p_name = "%s_prot.fasta" % region_name
    RegionSequencesFilename = "%s%s.fasta" % (
        Definitions.FastaRegionPrefix, region_name)

    RegionSequencesFilepath = os.path.join(options.WorkingDirectory,
                                           RegionSequencesFilename)

    if not os.path.isfile(RegionSequencesFilepath):
        print("Region sequences file not found at %s." %
              RegionSequencesFilepath)

    return list(SeqIO.parse(RegionSequencesFilepath, format="fasta"))


def EvaluateAllSequencesAllTranslationWindows(options,
                                              TemplateProtein,
                                              sequences):
    AllRegionSequences = []
    AllRecommendedWindows = []
    TotalSequences = 0
    for sequence in sorted(sequences, key=lambda s: s.id):
        if not EvaluateSequenceByUnknownBase(options, sequence):
            print("Discarding inaccurate sequence.")
            continue

        TotalSequences += 1
        (RecommendedWindow, WindowSequence, score) =\
            processAllTranslationWindows(
                options,
                TemplateProtein,
                sequence,
                options.Verbose
            )

        AllRecommendedWindows.append(RecommendedWindow)
        AllRegionSequences.append(WindowSequence)

        print("Correct Window: %i %i" % RecommendedWindow)
        print(">%s" % WindowSequence)
    return AllRecommendedWindows, AllRegionSequences, TotalSequences


def AnalyzeRegion(options, RegionSequenceSource):

    region_name = options.RegionName

    if not region_name:
        print("Region name undefined.")
        exit(1)

    sequences = LoadSequences(options, region_name)
    source_seq = RegionSequenceSource.fetchGeneSequence(region_name)

    if source_seq is None:
        return 0, 0

    TemplateProtein = SeqIO.SeqRecord(
        source_seq.translate(),
        id=region_name,
        description=""
    )

    AllRegionSequences = []
    for i in range(100):
        AllRecommendedWindows, AllRegionSequences, TotalSequences =\
            EvaluateAllSequencesAllTranslationWindows(options,
                                                      TemplateProtein,
                                                      sequences)

        if len(list(set(AllRecommendedWindows))) == 1:
            print("Found correct window.")
            break

    if options.WriteFiles:
        BuildOutputAlignments(
            options,
            region_name,
            AllRegionSequences,
            TemplateProtein
        )

    print("Rate for %s: %.2f%%" % (region_name, 0))
    return 0, 0


def BuildOutputAlignments(options, region_name,
                          AllRegionSequences, TemplateProtein):
    OutputProteinFilePrefix = os.path.join(options.WorkingDirectory,
                                           "Protein_%s" % region_name)

    OutputProteinFilePath = OutputProteinFilePrefix + ".fasta"
    OutputAlignmentFilePath = OutputProteinFilePrefix + ".aln"
    OutputTreeFilePath = OutputProteinFilePrefix + ".dnd"

    with open(OutputProteinFilePath, 'w') as f:
        SeqIO.write(AllRegionSequences, f, format="fasta")

    cmd = ClustalOmegaCommandline(Definitions.ClustalCommand,
                                  infile=OutputProteinFilePath,
                                  outfile=OutputAlignmentFilePath,
                                  guidetree_out=OutputTreeFilePath,
                                  outfmt="clustal",
                                  force=True)

    OutputProteinReferenceFilePath = os.path.join(
        options.WorkingDirectory,
        "Protein_ref_%s.fasta" % region_name
    )

    with open(OutputProteinReferenceFilePath, 'w') as f:
        SeqIO.write(TemplateProtein, f, format="fasta")

    cmd()


def runDirectory(options, RegionSequenceSource):
    WantedFileQuery = r"%s([\w\d]+).fasta" % Definitions.FastaRegionPrefix

    files = [
        f for f in os.listdir(options.WorkingDirectory)
        if re.findall(WantedFileQuery, f)
    ]

    if not files:
        print("No region fasta files detected.")
        exit(1)

    AllGapPresence = 0
    for f in files:
        region_name = re.findall(WantedFileQuery, f)[0]

        # Clone The Options (not quite a good development option)
        opt = copy.deepcopy(options)
        opt.RegionName = region_name

        successPercentage, HasGaps = AnalyzeRegion(opt, RegionSequenceSource)

        if options.WriteFiles:
            if successPercentage < 100:
                message = "%s with %.2f%%\n" % (region_name, successPercentage)
                with open("log", 'a') as f:
                    f.write(message)

        if HasGaps:
            AllGapPresence += 1

    print("Gaps presence: %.2f%%" % (100 * AllGapPresence / len(files)))


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", dest="RegionName")
    parser.add_argument("-d", "--dir", dest="WorkingDirectory", default=".")
    parser.add_argument("--norf", dest="OnlyOpenReadingFrame",
                        action="store_false", default=True)

    parser.add_argument("--discardimp",
                        dest="DiscardImperfectSequenceThreshold",
                        type=int, default=5)

    parser.add_argument("--nowrite",
                        dest="WriteFiles",
                        action="store_false",
                        default=True)

    parser.add_argument("-v", "--verbose", dest="Verbose", action="store_true")

    options = parser.parse_args()
    return options


def GetRegionSequenceSource(WorkingDirectory, genomeDirectory="genomes"):

    InformationFile = OutputFile.AnalysisInformation(WorkingDirectory)
    InformationFile.read()

    AnnotationPath = InformationFile.content["annotation"]

    GenomeFilePaths = genomeManager.readGenomeFolder(genomeDirectory)

    GenomeFeatures = list(SeqIO.parse(AnnotationPath, format="genbank"))

    return bfps.BruteForcePrimerSearcher(
        GenomeFeatures, GenomeFilePaths)


def Execute(options):
    RegionSequenceSource = GetRegionSequenceSource(options.WorkingDirectory)

    if options.RegionName:
        AnalyzeRegion(options, RegionSequenceSource)
    else:
        runDirectory(options, RegionSequenceSource)


def main():
    options = parse_arguments()
    Execute(options)


if __name__ == "__main__":
    main()
