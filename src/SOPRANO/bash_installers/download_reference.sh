#!/bin/bash

ref=$1
release=$2
data_dir=$3

if [ "$ref" == "GRCh37" ]
then
  fname="Homo_sapiens.GRCh37.dna.toplevel.fa"
  fname_gz="$fname.gz"
  link="https://ftp.ensembl.org/pub/grch37/release-$release/fasta/homo_sapiens/dna/$fname_gz"
elif [ "$ref" == "GRCh38" ]
then
  fname="Homo_sapiens.GRCh38.dna.toplevel.fa"
  fname_gz="$fname.gz"
  link="https://ftp.ensembl.org/pub/release-$release/fasta/homo_sapiens/dna/$fname_gz"
else
  echo "Genome ref $ref not configured"
  exit 1
fi

gz_path="$data_dir/$fname_gz"
fa_path="$data_dir/$fname"
fai_path="$fa_path.fai"
chrom_path="$data_dir/$(basename $fa_path .fa).chrom"

if [ ! -f "$fa_path" ]
then
    echo "Download $ref reference file? [Enter 'y' to download]"
    read response

    if [ "$response" == "y" ]
    then
        echo "downloading $link to $gz_path"
        wget "$link" -O "$gz_path" && echo "download complete."

        echo "decompressing $gz_path to $fa_path"
        gunzip -c "$gz_path" > "$fa_path" && echo "decompression complete."

        echo "removing compressed file $gz_path"
        rm "$gz_path"
    else
      echo "Interpreting $response as 'No'"
    fi
else
    echo "$ref reference already detected: $fa_path"
fi

if [ ! -f "$fai_path" ]
then
  echo "building fasta index file: $fai_path"
  samtools faidx "$fa_path" -o "$fai_path"
fi

if [ ! -f "$chrom_path" ]
then
  echo "computing chrom sizes: $chrom_path"
  cut -f1,2 "$fai_path" > "$chrom_path"
fi