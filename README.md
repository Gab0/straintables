## About

linkageMapper is a tool that helps evaluate the difference among gene loci across *Toxoplasma gondii* genomes.
<br>
It takes fully sequenced genomes for diverse strains of *T. gondii*, then find Amplicons inside a pair of Primers on those genomes.
The primers may be user-defined or found by the software when a locus name but no primers are defined.
<br>
Then it counts the SNPs that diverge from the amplicons found, and builds one Dissimilarity Matrix for each Amplicon/Gene/Primer_Pair that shows how the genomes are different at those regions.
<br>
More statistical analysis on the Dissimilarity Matrixes are carried, mostly using python's `skbio` module. The interpretation of analysis is under construction.
<br>
By looking at a pair of D. Matrixes at a time, both corresponding to locus that are neighbors, the user may have an idea of the recombination frequence of the studied organism.


## Setup

### Install Python Requirements

1. Clone this repository.

2. Install python dependencies: `$sudo pip install -r requirements.txt`

3. Make sure you have clustalW2 on your path.

### Fetch genomes and annotation files

The following code downloads all required genomes and annotations from NCBI,
populating the folders `genomes` and `annotations`. Standard usage requires this
execution only once.

```
$python linkageMapper/fetchDataNCBI.py
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


## Usage

1. Put the wanted Locus names, ForwardPrimers and ReversePrimers on a `.csv` file inside the `Primer` folder. The primer sequences are optional, leave blank to trigger the automatic primer search. Look for the examples.

2. `linkagePipeline.sh` is a bash script that organizes the workflow, calling main python scripts at `linkageMapper` module directory.

3. Check the results at the result folder that is equal to the `Primer` file selected for the run. Result folders are down the `Alignments` folder.


### Example 1

```
$ python linkageMapper/initializePrimerFile.py -i annotations/<chromosome annotation file name> -o Primers/TEST.csv
$ bash linkagePipeline.sh TEST
$ python linkageMapper/walkChromosomeResult.py -i Alignments/TEST
```

### Example 2

* First make your own primer file, with custom primers or blank primer fields. The primer file will be named `Primers/chr_X.csv` for this example.

```
$ bash linkagePipeline.sh chr_X
```

* Then view similarity matrixes and phylogenetic trees on `pdf` files at `Alignments/chr_X` folder.

### Example 3

* Make your own primer file, as explained in `Example 2`.
* The first locus should be one that is used commonly as *T. gondii* PCR-RFLP analysis
* For the primer files, we are using primers for SAG3 locus, and then other primers located at chromosome XII.

```
$ bash linkagePipelinea.sh chr_XII clustal SAG3
```

* Then view similarity matrixes and phylogenetic trees on `pdf` files at `Alignments/chr_X` folder.

## Results

* As the linkagePipeline.sh unfolds, `Alignments/[primer batch name]` folder will be created and populated with files of various kinds, in this order:

1. `.fasta` Sequence files, one holding the amplicon found for each loci.
2. `.aln` Alignment files, one for each loci.
3. `.aln.npy` Dissimilarity Matrix files, one for each loci.
4. `.pdf` Dissimilarity Matrix Plot files, one for each loci;


### Result Analysis Tools

Some python scripts on the main module are not called within linkagePipeline.sh. They are optional analysis tools and should be launched by the user.

1. `linkagePipeline/walkChromosomeResult.py` The basic one. This will build a slide presentation of plots, each plot showing a pair of similarity matrixes.
The matrix show will start with the first locus against the second locus chosen, in the order of the Primer file, along with some extra information. Then, the second locus will be compared to the third and so it goes on.



