#!/bin/bash

fasta_file="$DATA_DIR_PATH/ensemble_transcriptID.fasta"
zipped_fasta_file="$DATA_DIR_PATH/ensemble_transcriptID.fasta.gz"

if [ ! -f "$fasta_file" ]
then
  echo "-- decompressing transcript IDs"
  gunzip -k "$zipped_fasta_file"
  echo "-- decompressed: $fasta_file"
else
  echo "-- decompressed transcript IDs detected: $fasta_file"
fi