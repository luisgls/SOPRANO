#!/bin/bash

ref=$1
data_dir=$2

if [ "$ref" == "grch37" ]
then
  fname="Homo_sapiens.GRCh37.dna.toplevel.fa"
  link="https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
elif [ "$ref" == "grch38" ]
then
  fname=Homo_sapiens.GRCh38.dna.toplevel.fa
  link="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
else
  echo "Genome ref $ref not configured"
  exit 1
fi

gz_path="$data_dir/$fname.gz"
fa_path="$data_dir/$fname"

if [ ! -f "$fa_path" ]
then
    echo "Download $ref reference file? [Enter 'y' to download]"
    read response

    if [ "$response" == "y" ]
    then
        echo "downloading..."
        wget $link -o "$gz_path" && echo "download complete."

        echo "decompressing..."
        gunzip -c "$gz_path" > "$fa_path" && echo "decompression complete."

        rm "$gz_path"
        echo "$ref file written: $fa_path"
    else
      echo "Interpreting $response as 'No'"
    fi
else
    echo "$ref reference already detected: $fa_path"
fi
