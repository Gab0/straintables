# About

linkageMapper is a tool that helps evaluate the difference among gene loci across *Toxoplasma gondii* genomes.

It takes fully sequenced genomes for diverse strains of *T. gondii*, then find Amplicons inside a pair of primers on those genomes.

The primers may be user-defined or found by the software when a locus name is defined without attached primers.

Analysis proceeds while it counts the SNPs that diverge from the amplicons found at each locus, and builds one
 [dissimilarity matrix](https://en.wikipedia.org/wiki/Distance_matrix) for each Locus (which boundary is defined by a primer pair). 

The matrix can show how the genomes are different at those regions.

## Inside The Pipeline

### 1) Primer Docking, fetching the Amplicons

This step is carried on `linkageMapper/primerFinder.py`. For each designated loci, the app will try to find the complement and/or the original sequence
of both primers on all genomes. If both primers are found in a genome, the sequence between those primers is extracted and it proceeds to the next genome.
If every genome got its amplicon for the current locus, the locus was sucessfull and the script goes to the next one.

If for some reason not every genome is sucessfull with given pair of primers, the script downloads the locus sequence (if the locus name defined by the user matches a gene/locus on the genome annotation) from NCBI and randomizes one primer near the beginning of the gene sequence and other near the end.

Some available genomes are complement-reversed, and the app will make sure that loci sequences for all genomes are in the same orientation.


### 2) Amplicon Sequences Alignment

After getting the loci sequence from all the genomes, the visualization of the differences among genomes is done in two fronts:

#### 2a) Dissimilarity Matrix

1. The multifasta file containing sequence for one loci among all genomes is passed through ClustalW2
2. The the SNPs are detected and scored.
3. One Dissimilarity Matrix is created, showing which genome groups have similar locus.
4. Dissimilarity Matrixes can be viewed individually as `.pdf` files, `.npy` python files, or at the main visualization tool `walkChromosomeResult.py`.

#### 2b) MeshClust Clustering

1. The primary locus multifasta file is sent to MeshClust, which will detect clusters among genome's locus. Default MeshClust identity parameters is `0.999`.
2. The output of MeshClust is parsed at the visualization tool, which decorates genomes names at the Dissimilarity Matrix labels according to it's cluster group.

### 3) Visualization

Afther the pipeline executes the main scripts, the user can execute `linkageMapper/walkChromosomeResult.py` and load the folder in order to view the results.


![](https://raw.githubusercontent.com/Gab0/linkageMapper/master/walkChr.jpg?raw=true)

More statistical analysis on the Dissimilarity Matrixes are carried, mostly using python's `skbio` module. The interpretation of analysis is under construction.

By looking at a pair of D. Matrixes at a time, both corresponding to locus that are neighbors, the user may have an insight on data of the studied organism, like the recombination frequency.

# Setup

### Install Python Requirements

1. Clone this repository.

2. Install python dependencies: `$sudo pip install -r requirements.txt`

### Fetch genomes and annotation files

The following code downloads all required genomes and annotations from NCBI,
populating the folders `genomes` and `annotations`. Standard usage requires one time execution of the following command:

```
To download defaults: Toxoplasma gondii genomes & TGME49 annotations;
$python linkageMapper/fetchDataNCBI.py 

Same result as above:
$python linkageMapper/fetchDataNCBI.py --organism "Toxoplasma gondii" --strain ME49

With lactobacillus plantarum & WCFS1 annotations:
$python linkageMapper/fetchDataNCBI.py --organism "Lactobacillus plantarum" --strain WCFS1
```

Instead, you can manually add desired genomes and annotations, as explained in the next subsections:

#### Annotations

* The annotation file serves as a guide for the analysis, as they contain most loci names for *T. gondii*.
* Recommended annotation is the one for TGME49, which at the moment seems to be the most complete one.
* All `.gff` annotation files, one for each chromosome, should be placed at `annotations` folder.

#### Genomes

* The genome files are the root of the analysis.
* One multifasta file per strain.
* They should be placed at the `genomes` folder.


## Required Software

### ClustalW2

The alignment step of `linkageMapper` requires [ClustalW2](http://www.clustal.org/clustal2/) installed on your
system.

### MeShClust [optional/redundant]

The recombination analysis step of `linkageMapper` has [MeShCluSt](https://github.com/TulsaBioinformaticsToolsmith/MeShClust) as an optional dependency.

Having it installed on the system will enable genome group clustering to be totally independend from the alignment software, as MeShCluSt does the clustering
on top of unaligned `.fasta` files.

By not having it, the clustering will be made off distance matrix information only. 

At this point, it is not important to have MeShCluSt installed.


# Usage

1. Put the wanted Locus names, ForwardPrimers and ReversePrimers on a `.csv` file inside the `Primer` folder. The primer sequences are optional, leave blank to trigger the automatic primer search. Look for the examples.

2. `linkagePipeline.sh` is a bash script that organizes the workflow, calling main python scripts at `linkageMapper` module directory.

3. Check the results at the result folder that is equal to the `Primer` file selected for the run. Result folders are down the `Alignments` folder.


#### Example 1: Automatic Locus Selection with Automatic Primer Search.

```
$ python linkageMapper/initializePrimerFile.py -i annotations -c X -o Primers/TEST.csv -p 0.01
$ bash linkagePipeline.sh TEST
$ python linkageMapper/walkChromosomeResult.py -i Alignments/TEST
```



#### Example 2: Custom Locus Selection with Automatic Primer Search

* Make your own primer `.csv` file, named `Primers/chr_X.csv` for this example. It should have blank primer fields. 

```
<@file: Primers/chr_X.csv>
LocusName,ForwardPrimer,ReversePrimer
CDPK,,
IMC2A,,
AP2X1,,
TGME49_227830,,
```


Then, execute:
```
$ bash linkagePipeline.sh chr_X
```

* Then view similarity matrixes and phylogenetic trees on `pdf` files at `Alignments/chr_X` folder.


#### Example 3: Custom Loci Selection, Custom Primer Search

* Follow Example 2, except now the primer file can have a pair of primers designed for each loci:
* Some primers, if missing or problematic, will trigger the automatic primer search.

```
<@file: Primers/chr_X.csv>
LocusName,ForwardPrimer,ReversePrimer
CDPK1,ACAAAGGCTACTTCTACCTC,TTCTATGTGGGGATGCAGAG
IMC2A,,GACGGACGCATGGCTTGCTG
AP2X1,GCTCAAGCTGCTCCCCGGGC,TCGACGGAGGTGCTCCAACC

```

<!--- DEPRECATED?
### Example 4:

* Make your own primer file, as explained in `Example 2`.
* The first locus should be one that is used commonly as *T. gondii* PCR-RFLP analysis
* For the primer files, we are using primers for SAG3 locus, and then other primers located at chromosome XII.

```
$ bash linkagePipelinea.sh chr_XII clustal SAG3
```

* Then view similarity matrixes and phylogenetic trees on `pdf` files at `Alignments/chr_X` folder.
-->

# Results

* As the linkagePipeline.sh unfolds, `Alignments/[primer batch name]` folder will be created and populated with files of various kinds, in this order:

1. `.fasta` Sequence files, one holding the amplicon found for each loci.
2. `.aln` Alignment files, one for each loci.
3. `.aln.npy` Dissimilarity Matrix files, one for each loci.
4. `.pdf` Dissimilarity Matrix Plot files, one for each loci;


## Result Analysis Tools

Some python scripts on the main module are not called within linkagePipeline.sh. They are optional analysis tools and should be launched by the user.

1. `linkagePipeline/walkChromosomeResult.py` The basic one. This will build a slide presentation of plots, each plot showing a pair of similarity matrixes.
The matrix show will start with the first locus against the second locus chosen, in the order of the Primer file, along with some extra information. Then, the second locus will be compared to the third and so it goes on.



