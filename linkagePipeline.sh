#!/bin/bash


PRIMER_CODE=$1

OUTPUT_DIR="Alignments/${PRIMER_CODE}"

mkdir -p "${OUTPUT_DIR}"

if python linkageMapper/primerFinder.py -r Primers/"${PRIMER_CODE}".csv -o "${OUTPUT_DIR}" -p;
then
    echo Proceeding analysis...
else
    echo Failure.
    exit
fi

# THIS GETS ALL LOCI NAMES FROM .csv PRIMER DESCRIPTOR!
# WORKS LIKE THIS -> LOCI=(SAG2 ROP32 BIN3)
readarray -t LOCI < <(cat Primers/"${PRIMER_CODE}".csv|awk -F"," '{print $1}'|tail -n +2)


ALNMODE=$2

CLONAL=$3

if [ -z $CLONAL ]
then
    EXPLICT_CLONAL=""
else
    EXPLICIT_CLONAL="-s 'vs clonal ${CLONAL}'"
fi

for _L in "${LOCI[@]}"
do
    OUTPUT_FILE_PREFIX="${OUTPUT_DIR}/LOCI_${_L}"
    python linkageMapper/makeMultifasta.py -l "${_L}" -r "${LOCI[0]}" -i "${OUTPUT_DIR}"/Sequences.csv

    if [ $ALNMODE = "tcoffee" ];
    then
        tcoffee -in "${OUTPUT_FILE_PREFIX}".fasta -out "${OUTPUT_DIR}"
    fi

    if [ $ALNMODE = "clustal" || -z $ALNMODE];
    then
        clustalw2 -INFILE="${OUTPUT_FILE_PREFIX}".fasta -OUTFILE="${OUTPUT_FILE_PREFIX}".aln
    fi

    if [ $ALNMODE = "muscle" ];
    then
        muscle -in "${OUTPUT_FILE_PREFIX}".fasta -out "${OUTPUT_FILE_PREFIX}".aln
    fi

    # BUILD TREE;
    clustalw2 -INFILE="${OUTPUT_FILE_PREFIX}".aln -tree

    # POST-ALIGN ANALYSIS;
    python linkageMapper/drawTree.py "${OUTPUT_FILE_PREFIX}".ph
    python linkageMapper/detectMutations.py -i "${OUTPUT_FILE_PREFIX}".aln "${EXPLICIT_CLONAL}"
done

#SIMILARITY MATRIX DIFFERENCES;
python linkageMapper/compareHeatmap.py -d "${OUTPUT_DIR}"

# ANALYZE MATRIX DIFFERENCES
python linkageMapper/matrixAnalysis.py -d "${OUTPUT_DIR}"


# CLEANUP USELESS FILES;
rm "${OUTPUT_DIR}"/*.dnd
# rm "${OUTPUT_DIR}"/*.aln.npy
rm "${OUTPUT_DIR}"/*.ph
rm "${OUTPUT_DIR}"/*.aln.csv
