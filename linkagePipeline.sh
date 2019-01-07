#!/bin/bash

# SHOW BEAUTIFUL ASCII ART;
cat logo.txt

PRIMER_CODE=$1

# DEFINE OUTPUT DIRECTORY;
OUTPUT_DIR="analysisResults/${PRIMER_CODE}"

# MAKE SURE OUTPUT DIRECTORY EXISTS;
mkdir -p "${OUTPUT_DIR}"

# INITIALIZE DEFAULT COMMAND LINE OPTIONS;
DO_AMPLICON=1
DO_ALIGNMENT=1
ALNMODE="clustal"
DO_MESHCLUST=1

# FETCH COMMAND LINE ARGUMENTS;
ALL_ARGS=($1 $2 $3 $4 $5 $6 $7)

# PARSE COMMAND LINE ARGUMENTS;
for OPT in "${ALL_ARGS[@]}"
do

    if [ $OPT = "noamplicon" ];
    then
        DO_AMPLICON=0
    fi

    if [ $OPT = "noalign" ];
    then
        DO_ALIGNMENT=0
    else
        for _ALNMODE in "${ALL_ALNMODES[@]}"
        do
            if [ $_ALNMODE = $OPT ];
            then
                ALNMODE=$OPT
            fi
        done
    fi
done



# OUTPUT INFO;
echo "Alignment mode is $ALNMODE"
echo ""

if meshclust;
then
    echo "MeShClust enabled!"
else
    echo "MeShClust not found!"
    DO_MESHCLUST=0
fi

echo '\n'
echo "Running pipeline..."
echo ""

# RUN AMPLICON SEARCH IF NOT DISABLED BY OPTIONS;
if [ $DO_AMPLICON = 1 ];
then
    if python linkageMapper/primerFinder.py \
              -i Primers/"${PRIMER_CODE}".csv \
              -o "${OUTPUT_DIR}" -p \
              -r "${LOCI[0]}"
    then
        echo "Proceeding analysis..."
    else
        echo "Failure."
        exit 1
    fi
else
    echo "No amplicon Mode!"
fi

# THIS GETS ALL LOCI NAMES FROM .csv PRIMER DESCRIPTOR!
# WORKS LIKE THIS -> LOCI=(SAG2 ROP32 BIN3)
readarray -t LOCI < <(cat "${OUTPUT_DIR}"/MatchedPrimers.csv|awk -F"," '{print $1}'|tail -n +2)

CLONAL=$3

if [ -z $CLONAL ];
then
    EXPLICT_CLONAL=""
else
    EXPLICIT_CLONAL="-s 'vs clonal ${CLONAL}'"
fi

# ITERATE ALL LOCI FROM MatchedPrimers.csv;
for _L in "${LOCI[@]}"
do
    OUTPUT_FILE_PREFIX="${OUTPUT_DIR}/LOCI_${_L}"


    # RUN ALIGNMENT IF NOT DISABLED BY OPTIONS;
    echo "Running alignment for ${OUTPUT_FILE_PREFIX}"
    echo ""
    if [ $DO_ALIGNMENT = 1 ];
    then
        if [ $ALNMODE = "tcoffee" ];
        then
            tcoffee -in "${OUTPUT_FILE_PREFIX}".fasta -out "${OUTPUT_DIR}"
        fi

        if [ $ALNMODE = "clustal" ] || [ -z $ALNMODE ];
        then
            if timeout 60 clustalw2 \
                -INFILE="./${OUTPUT_FILE_PREFIX}".fasta \
                -OUTFILE="./${OUTPUT_FILE_PREFIX}".aln;
            then
                echo "Alignment Done!"
            else
                echo "Alignment Fails!"
                continue
            fi

            #2>>"${OUTPUT_DIR}/clustal_warnings.txt"
        fi

        if [ $ALNMODE = "muscle" ];
        then
            muscle -in "${OUTPUT_FILE_PREFIX}".fasta -out "${OUTPUT_FILE_PREFIX}".aln
        fi
    fi

    # BUILD TREE;
    echo "Building tree...."
    if clustalw2 -INFILE="./${OUTPUT_FILE_PREFIX}".aln -tree;
    then
        echo "Tree built."
        echo ""
        echo ""
    else
        echo "Tree fails."
        exit 1 
    fi

    # DRAW PHYLOGENETIC TREES FROM .ph CLUSTALW2 FILES;
    python linkageMapper/DrawGraphics/drawTree.py \
           -i "${OUTPUT_FILE_PREFIX}".ph \
           -o "${OUTPUT_FILE_PREFIX}".pdf

    # INTERPRETE MUTATIONS AS DISSIMILARITY MATRIXES;
    python linkageMapper/detectMutations.py \
           -i "${OUTPUT_FILE_PREFIX}".aln \
           "${EXPLICIT_CLONAL}"

    # MESHCLUST CLUSTERS;
    if [ $DO_MESHCLUST == 1 ];
       then
           meshclust "${OUTPUT_FILE_PREFIX}".fasta \
                     --output "${OUTPUT_FILE_PREFIX}".clst \
                     --id 0.999 --align
    fi
done

# SIMILARITY MATRIX DIFFERENCES;
python linkageMapper/compareHeatmap.py -d "${OUTPUT_DIR}"

# ANALYZE MATRIX DIFFERENCES
python linkageMapper/matrixAnalysis.py -d "${OUTPUT_DIR}"


# CLEANUP USELESS FILES;
rm "${OUTPUT_DIR}"/*.dnd
# rm "${OUTPUT_DIR}"/*.aln.npy
rm "${OUTPUT_DIR}"/*.ph
rm "${OUTPUT_DIR}"/*.aln.csv
