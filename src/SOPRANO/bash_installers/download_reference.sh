#!/bin/bash

ref=$1
release=$2
data_dir=$3

if [ "$ref" == "GRCh37" ]
then
  baseurl="https://ftp.ensembl.org/pub/grch37/release-$release/fasta/homo_sapiens/dna"
elif [ "$ref" == "GRCh38" ]
then
  baseurl="https://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna"
else
  echo "Genome ref $ref not configured"
  exit 1
fi

# Define filenames, urls and target paths

# Basenames
toplevel_basename="Homo_sapiens.$ref.dna.toplevel"
primary_basename="Homo_sapiens.$ref.dna.primary_assembly"

# fasta files
toplevel_fa="$toplevel_basename.fa"
toplevel_fa_path="$data_dir/$toplevel_fa"
primary_fa="$primary_basename.fa"
primary_fa_path="$data_dir/$primary_fa"

# compressed fasta files
toplevel_gz="$toplevel_fa.gz"
toplevel_gz_path="$data_dir/$toplevel_gz"
primary_gz="$primary_fa.gz"
primary_gz_path="$data_dir/$primary_gz"

# fasta index files
toplevel_fai="$toplevel_fa.fai"
toplevel_fai_path="$data_dir/$toplevel_fai"

# chrom files
toplevel_chrom="$toplevel_basename.chrom"
toplevel_chrom_path="$data_dir/$toplevel_chrom"

# source urls
toplevel_url="$baseurl/$toplevel_gz"
primary_url="$baseurl/$primary_gz"


# Download compressed toplevel fasta file if decompressed version not found
if [ ! -f "$toplevel_fa_path" ]
then
    echo "Download $ref toplevel reference file? [Enter 'y' to download]"
    read response

    if [ "$response" == "y" ]
    then
        echo "downloading $toplevel_url to $toplevel_gz_path"
        wget "$toplevel_url" -O "$toplevel_gz_path" && echo "download complete."

        echo "decompressing $toplevel_gz_path to $toplevel_fa_path"
        gunzip -c "$toplevel_gz_path" > "$toplevel_fa_path" && echo "decompression complete."

        echo "removing compressed file $toplevel_gz_path"
        rm "$toplevel_gz_path"
    else
      echo "Interpreting $response as 'No'"
    fi
else
    echo "$ref toplevel reference already detected: $toplevel_fa_path"
fi

# Download compressed primary fasta file if decompressed version not found
if [ ! -f "$primary_fa_path" ]
then
    echo "Download $ref primary assembly reference file? [Enter 'y' to download]"
    read response

    if [ "$response" == "y" ]
    then
        echo "downloading $primary_url to $primary_gz_path"
        wget "$primary_url" -O "$primary_gz_path" && echo "download complete."

        echo "decompressing $primary_gz_path to $primary_fa_path"
        gunzip -c "$primary_gz_path" > "$primary_fa_path" && echo "decompression complete."

        echo "removing compressed file $primary_gz_path"
        rm "$primary_gz_path"
    else
      echo "Interpreting $response as 'No'"
    fi
else
    echo "$ref primary assembly reference already detected: $primary_fa_path"
fi

# Run samtools faidx against fasta file if fasta index not found
if [ ! -f "$toplevel_fai_path" ]
then
  echo "building fasta index file: $toplevel_fai_path"
  samtools faidx "$toplevel_fa_path" -o "$toplevel_fai_path"
else
    echo "fasta index file already detected: $toplevel_fai_path"
fi

# Compute chrom sizes from fasta index file if not found
if [ ! -f "$toplevel_chrom_path" ]
then
  echo "computing chrom sizes: $toplevel_chrom_path"
  cut -f1,2 "$toplevel_fai_path" > "$toplevel_chrom_path"
else
    echo "chrom sizes file already detected: $toplevel_chrom_path"
fi