# SOPRANO: Selection On PRotein ANnotated regiOns
SOPRANO was developed to analyse selection in specific regions of the genome. It uses a VEP annotated file to estimate ON target dN/dS values and OFF target dN/dS values.

## Installation

#### To install first create a directory called SOPRANO and then clone the tool to this directory

```{bash}
mkdir /my/home/directory/SOPRANO
cd /my/home/directory/SOPRANO

git clone https://github.com/luisgls/SOPRANO.git .
```
#### Edit the head of the master script: run_localSSBselection_v3.sh
- Specify the basedirectory of the installation
BASEDIR=/my/home/directory/SOPRANO/

- Copy or link the genome file (e.g. hg19.genome) and the fasta file (e.g. hg19.fasta) to /my/home/directory/SOPRANO/data/

#### Now you should be able to run the tool if all the dependencies are met.

### Dependencies
- bedtools 2.26.0
- R-3.3.3 or higher.
- R library tidyr
- perl 5
- GNU command line tools
- Ensembl variant effect predictor v89 or higher (VEP)

#### Important Notes
- earlier versions of bedtools will not work
- tab encoding should be \t (might be a problem for windows/OSX versions)
- genome file is a two column file specifying the fasta id and the length of the sequence (see how to obtain it at the bottom)
- Restrict your input dataset to chromosomes 1-22 and X and Y. Remove the rest.

## Input file
The input file is the standard output of variant effect predictor using the following command line (by providing to vep the ensembl default input file format)
perl variant_effect_predictor.pl -i input -o input.annotated --cache --all_refseq --assembly GRCh37 --pick --symbol --no_stats --fasta genome.fasta
If you want to filter putative germline variants use the option --plugin ExAC when running VEP.

Example input files can be found on synapse: ID syn11681983

#### Important points before running
  - a) No header needed for input VEP file
  - b) VEP annotated first column must be in the format (chr_pos_ref/alt)
  - c) VEP annotated file must only have chromosomes that are 1,2,3,4...22 or uppercase X,Y
  - d) After you run vep with the option for ExAC frequencies, it would be necessary to remove all variants present in more than 0.1 percent of the population. You could apply the   filter using:
filter_vep -i input.annotated -f "ExAC_AF < 0.1 or not ExAC_AF" --ontology --filter "Consequence is coding_sequence_variant" 
  - e) Be sure that you are using the GNU command line if you are running in a MacOS (https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/)
  - f) Add dependencies to your path for easy running or hardcode the scripts
  - g) The genome file used in SSB is a two column file that contains the info of the name of the fasta id (column 1) and the length of that sequence (column 2).
  - h) The UNIX system used should be able to recognize \t as a tab separator

## Genomes
To get hg19 fasta genome, you can download it from UCSC:

```{bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

```

