## About

linkageMapper is a tool that helps evaluate the difference among gene loci across *Toxoplasma gondii* genomes.

## Setup

### Requirements
1. Clone this repository.

2. Install python dependencies: `$sudo pip install -r requirements.txt`

3. Make sure you have clustalW2 on your path.

### Fetch annotation files

* The annotation file serves as a guide for the analysis, as they contain most loci names for *T. gondii*.
* Recommended annotation is the one for TGME49, which at the moment seems to be the most complete one.
* Annotation file names should be like `genomes/chromosome_X` and onwards.

```
place at annotations folder.
TBD.

```

### Fetch genomes

* The genomes file are the root of the analysis.
* They should be placed at the `genomes` folder.

```

TBD.

```


## Usage

1. Put the wanted Locus names, ForwardPrimers and ReversePrimers on a `.csv` file inside the `Primer` folder. The primer sequences are optional, leave blank to trigger the automatic primer search. Look for the examples.

2. `linkagePipeline.sh` is a bash script that organizes the workflow, calling main python scripts at `linkageMapper` module directory.

3. Check the results at the result folder that is equal to the `Primer` file selected for the run. Result folders are down the `Alignments` folder.


### Example 1

```

$ python linkageMapper/initializePrimerFile.py -i annotations/chromosome_X -o Primers/walk_chromosome_X
$ bash linkagePipeline.sh walk_chromosome_X
$ python linkageMapper/walkChromosomeResult.py -i Alignments/walk_chromosome_X

```
### Example 2

* First make your own primer file, with custom primers or blank primer fields. That's the file `Primers/chr_X` for this example.

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

TBD.
