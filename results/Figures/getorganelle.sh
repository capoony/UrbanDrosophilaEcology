#!/bin/sh

#PBS -N getorganelle_surf
#PBS -o /media/inter/geuth/Garra/genome_assemblies/log/getorganelle_surf.log
#PBS -j oe
#PBS -l select=1:ncpus=35:mem=100g

# ==============================================================
# Extract mtDNA from short reads using GetOrganelle
# Annotate the fasta mtDNA using custom Python script
# Use the gb file to produce a PNG image of the annotated mtDNA
# ==============================================================

# ----------------------------------------------------------------------------------------------------------------------------------
# DATA PREPARATION, SETTING OF DIRECTORIES & MODULES LOADING

# Parse data and species variable from the config file
CONFIG="/media/inter/geuth/Garra/genome_assemblies/scripts/config_files/getorganelle.txt"

# Set Species variable

SPECIES=$( grep "^SPECIES=" "$CONFIG" | cut -d "=" -f2 )

if [ -n "$SPECIES" ]; then
    echo "Species: $SPECIES"
else
    echo "Species variable was not set correctly. Exiting"
    exit 1
fi

# Parse the DNBseq (R1,R2) files
DNB1=$( grep "^DNB1_${SPECIES}=" "$CONFIG" | cut -d "=" -f2 )
DNB2=$( grep "^DNB2_${SPECIES}=" "$CONFIG" | cut -d "=" -f2 )

if [ -n "$DNB1" ] && [ -n "$DNB2" ]; then
    echo "Using R1: $DNB1 and R2: $DNB2"
else
    echo "At least one of the two R1 and R2 files for $SPECIES were not parsed correctly. Exiting"
    exit 1
fi

# Parse the reference genome file and database
REFERENCE=$( grep "^REFERENCE_${SPECIES}=" "$CONFIG" | cut -d "=" -f2 )

# Load dependencies
source /media/inter/geuth/.venv/bin/activate
source /opt/anaconda3/etc/profile.d/conda.sh
conda activate getorganelle 

# Define output directories
OUT_DIR=/media/inter/geuth/Garra/genome_assemblies/results/${SPECIES}/assembly/mitochondrion/getorganelle
PYSCRIPT=/media/inter/geuth/Garra/genome_assemblies/scripts/pyscripts/mtgenome_annotation.py

mkdir -p "$OUT_DIR"
# ----------------------------------------------------------------------------------------------------------------------------------

# 1. Get mt reads from the short read R1 and R2 files
get_organelle_from_reads.py -1 "$DNB1" -2 "$DNB2" \
    --max-reads 3E8 -k 21,31,55,85,115 -F animal_mt \
    -t 125 \
    --continue \
    -o "$OUT_DIR" 

if [ ! -f "$OUT_DIR/${DB}.K115.complete.graph1.1.path_sequence.fasta" ]; then
    echo "Something went wrong. Analysis failed. Exiting"
    exit 1
else
    echo "Analysis completed successfully. Output files are in $OUT_DIR"
fi

mv "$OUT_DIR/${DB}.K115.complete.graph1.1.path_sequence.fasta" "$OUT_DIR/${SPECIES}_mtDNA.fasta"

# 2. Annotate the mtDNA with biopython
python3 "$PYSCRIPT" -f "$REFERENCE" -i "$OUT_DIR/${SPECIES}_mtDNA.fasta" -d "$SPECIES" -o "$OUT_DIR/${SPECIES}_mtDNA.gb"

if [ ! -f "$OUT_DIR/${SPECIES}_mtDNA.gb" ]; then
    echo "Annotation failed. Exiting"
    exit 1
else
    echo "Annotation completed successfully. Output file is $OUT_DIR/${SPECIES}_mtDNA.gb"
fi 