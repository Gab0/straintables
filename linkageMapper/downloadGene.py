#!/bin/python

from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
Entrez.email = "thegondii@mail.com"


def downloadFasta(geneName):
    Query = geneName + "[Gene] Toxoplasma gondii[Organism]"
    h = Entrez.esearch(db='gene', term=Query)
    h = Entrez.read(h)
    print(h.keys())
    ID = h['IdList'][0]

    print(ID)

    h = Entrez.efetch(db='gene', id=ID, retmode='text', rettype='fasta').read()
    h = Entrez.efetch(db='gene', id=ID, retmode='text', rettype='fasta').read()
    print(h)
    #h = Entrez.read(h)


def retrieveGeneSequence(genomeFeatures, geneName):
    for g, FeatureGroup in enumerate(genomeFeatures):
        for feature in FeatureGroup.features:
            A = FeatureGroup
            if feature.type == "gene":
                MATCH = False
                if "gene" in feature.qualifiers.keys():
                    if geneName in feature.qualifiers['gene']:
                        MATCH = True
                else:
                    if geneName in feature.qualifiers['locus_tag']:
                        MATCH = True
                if MATCH:
                    return FeatureGroup.description, feature.location

def fetchSequence(location, chr_descriptor, genomeFilePath):
    genome = list(SeqIO.parse(genomeFilePath, format="fasta"))

    for c, Chromosome in enumerate(genome):
        if chr_descriptor in Chromosome.description:
            Sequence = Chromosome.seq[location.start.position:location.end.position]
            if location.strand == -1:
                Sequence = Sequence.reverse_complement()
            return Sequence
        # print(dir(Chromosome))
        # print(Chromosome.name)
        # print(Chromosome.id)
        # print(Chromosome.description)
