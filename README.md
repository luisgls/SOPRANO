# SOPRANO: Selection On PRotein ANnotated regiOns
SOPRANO method was developed to quantify selection in specific regions of the 
genome (Zapata et al, in revision). It calculates ON- and OFF-target dN/dS 
using a set of annotated somatic point mutations and a transcript coordinates 
file. 

- The set of mutations (missense/truncating) and their respective functional 
annotation using ensemblVEP. Mutations can be filtered a priori by the user 
(i.e. by predicted binding affinity or by expression status).

- The set of ensembl transcript coordinates where selection will be estimated. 
ON-dN/dS is the value calculated inside the coordinates provided using a 
192-trinucleotide correction signature obtained "on-the-fly" from the input 
mutation file. Alternatively, the user can provide a pre-calculated 
trinucleotide mutation frequency file. Importantly, ON dN/dS and OFF dN/dS 
(the portion outside the coordinates provided) will be calculated only using 
transcripts defined in this file. 

## Installation

SOPRANO currently requires Mac or Linux OS. All dependencies (including 
bioinformatics tools such as VEP and bedtools) are built within a conda 
environment for coherent reproducibility. The suite of SOPRANO tools is 
assembled as a Python package, which is installed into the wider conda 
environment.

To install the environment and SOPRANO code, users are therefore required 
to have a conda (or mamba) installation available. From the repository root
users should run 

```. setup.sh```

which will construct and activate the `soprano-dev` environment. After the 
first installation, users can activate the pre-built environment with

``` conda activate soprano-dev ```

## Running SOPRANO

Once the `soprano-dev` environment has been installed and activated, the 
following options are available.

### [Option 1] Command line interface

SOPRANO can be run from the command line via:

``` RUN_SOPRANO <arguments> ```

**Required arguments**

Definitions should follow the short (`-s`) or long (`--long`) argument flags.

- `-n | --name`

    The name of the pipeline run.

- `-o | --output`

    The path to the output file directory (i.e. the cache) for the SOPRANO
    pipeline run.

- `-i | --input`

    The path to the input VEP annotated file that defines mutations.

- `-b | --bed_file`

    The path to the bed file containing protein coordinates named by their 
    transcript ID (e.g. ENSTXXXXXX 123 135)

**Optional arguments**

Definitions should follow the short (`-s`) or long (`--long`) argument flags.

- `-t | --target_regions`

    The path to the bed file that defines regions for randomization.
    Defaults to None - implying 

- `-s | --seed`

    The seed value to be used in randomization process. If defined, users 
    should be able to reproduce their results when randomization is applied.
    This value should be an integer greater than zero. Defaults to "-1", 
    implying no seed value chosen.

- `-r | --reference`

    The name of the genome reference. Currently this CLI supports GRCh37 and 
    GRCh38, which may be passed as arguments. By default, GRCh37 is chosen.

- `-q | --release`

    This corresponds to the Ensembl release number of the genome reference
    file. By default, the CLI will select release 110.

**Optional switches**

These switches are activated by passing the long (`--long`) flag to the command
line with no subsequent argument. They alter the runtime behaviour of SOPRANO.

- `--use_random`

    Run SOPRANO with randomization methods.

- `--use_ssb192`

    Run SOPRANO using SSB192 substitutions, otherwise will use SSB7.

- `--exclude_drivers`

    Run SOPRANO with driver genes excluded. See drivers at 
    `/src/SOPRANO/data/genes2exclude.txt`.

### [Option 2] Local web interface

This feature is in development (fall 2023). Users will be able to run SOPRANO 
from their web browser without the command line interface.

## Downloading the genome reference files

Once the `soprano-dev` environment has been installed and activated, users can
download the necessary genome reference files for their analyses via 

```GET_GENOMES -r <reference>```

Currently supported references are `GRCh37` and `GRCh38`. By default, the 
download willl correspond to the ensembl 110 data release(s). Users can 
additionally prescribe a different release if they wish, e.g.

```GET_GENOMES -r GRCh37 --release 109```

Users will be prompted to download the toplevel and primary assembly genome
files. The former is strictly required to run SOPRANO, whilst the latter will
accelerate the annotation procedure (currently in development).

Once downloaded, SOPRANO will automatically compute the corresponding fasta
index files and chromosome sizes.

The downloaded data mimics the cache structure of ensembl VEP, i.e., files will
be downloaded to `/src/SOPRANO/data/homo_sapiens/<release>_<reference>/`.

**Warning:** Downloading all the files associated with a given reference will
consume ~70GB of memory. Users are therefore advised to only download the 
primary assembly files if they are required for annotations.

[//]: # (#### Create a directory and clone the repository)

[//]: # ()
[//]: # (```{bash})

[//]: # (mkdir tools)

[//]: # (cd tools)

[//]: # ()
[//]: # (git clone https://github.com/luisgls/SOPRANO.git )

[//]: # (```)

[//]: # (### Install dependencies)

[//]: # (- &#40;For MacOSX users&#41; GNU command line tools &#40;https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/&#41;)

[//]: # (- bedtools 2.26.0 or later &#40;https://bedtools.readthedocs.io/en/latest/content/installation.html&#41;)

[//]: # (- R-3.3.3 or later &#40;https://www.r-project.org/&#41;. For quick install if you have brew &#40;https://formulae.brew.sh/formula/r&#41;)

[//]: # (- R library tidyr &#40;https://tidyr.tidyverse.org/&#41;)

[//]: # (- perl 5 &#40;https://www.perl.org/get.html&#41;)

[//]: # (- Ensembl variant effect predictor v89 or higher &#40;VEP&#41; &#40;https://www.ensembl.org/info/docs/tools/vep/index.html&#41;)

[//]: # ()
[//]: # (#### Copy a local version and edit the master script to run locally &#40;run_localSSBselection_v4.sh&#41;)

[//]: # ()
[//]: # (```{bash})

[//]: # (cp run_localSSBselection_v4.sh run_localSSBselection_v4_local.sh)

[//]: # (```)

[//]: # ()
[//]: # (- Specify the base directory of the installation)

[//]: # (BASEDIR=/the/home/directory/SOPRANO/)

[//]: # ()
[//]: # (- Specify a tmp directory)

[//]: # (TMP=/my/tmp/)

[//]: # ()
[//]: # (- Specify the reference genome sequence fasta file &#40;e.g. hg19.fasta&#41;)

[//]: # (FASTA=/my/directory/to/hg19.fasta)

[//]: # ()
[//]: # (- Specify the reference genome length file &#40;e.g. hg19.genome&#41;)

[//]: # (GENOME=/my/directory/to/hg19.genome)

## Input mutation files 

This data corresponds to the SOPRANO CLI argument `-i | --input`

- Mutations should be annotated with VEP using the default output format 
(https://www.ensembl.org/info/docs/tools/vep/vep_formats.html).
The input file is the standard output of variant effect predictor using the 
following command line (by providing to vep the ensembl default input file 
format) ```perl variant_effect_predictor.pl -i input -o input.annotated --cache --all_refseq --assembly GRCh37 --pick --symbol --no_stats --fasta genome.fasta```
If you want to filter putative germline variants use the option --plugin ExAC
when running VEP. It is important that you restrict your analysis to the list
of ensemb transcripts observed in data folder, be aware of updates on the 
EnsemblID from previous versions.

- VCF file annotated with VEP using VCF format
To convert a VCF annotated file to the input format for SOPRANO, the user can
run vep-annotation-reporter from 
(https://vatools.readthedocs.io/en/latest/vep_annotation_reporter.html) using
the following command:

```{bash}
vep-annotation-reporter -o OUTPUT.tsv INPUT.vcf.annotated Allele Gene Feature Feature_type Consequence cDNA_position CDS_position Protein_position Amino_acids Codons Existing_variation 
```
Example input files can be found on synapse: ID syn11681983

## Input Immunopeptidome files

This data corresponds to the SOPRANO CLI argument `-b | --bed_file`

- A bed file containing the transcript name (ensembl transcript identifier), 
  the start and end of the peptide that binds the MHC complex

- Example immunopeptidome files can be found under the immunopeptidome
  directory.

## Output file

The output of SOPRANO is a table report, with headers

```{bash}
coverage ON_dnds ON_lowci ON_highci ON_muts OFF_dnds OFF_lowci OFF_highci OFF_muts P-val ON_na ON_NA ON_ns ON_NS OFF_na OFF_NA OFF_ns OFF_NS
```

- ```coverage``` Two options: ExonicOnly and ExonicIntronic. The latest should be used if there are intronic mutations in the mutation file. The algorithm uses intronic mutations to improve the background counts of silent mutations.

- ```ON_dnds``` dN/dS of the target region provided in the bed file

- ```ON_lowci``` lower value for the 95% CI of the target

- ```ON_highci``` upper value for the 95% CI of the target

- ```ON_muts``` number of mutations observed inside the target region

- ```OFF_dnds``` dN/dS of the OFF-target region provided in the bed file

- ```OFF_lowci``` lower value for the 95% CI of the OFF-target

- ```OFF_highci``` upper value for the 95% CI of the OFF-target

- ```OFF_muts``` number of mutations observed outside the target region

- ```P-val``` P-value estimated from the comparison of the confidence intervals from ON and OFF dN/dS values

- ```ON_na``` Observed number of nonsilent mutations ON target

- ```ON_NA``` Number of nonsilent sites (corrected) ON target

- ```ON_ns``` Observed number of silent mutations ON target

- ```ON_NS``` Number of silent sites (corrected) ON target

- ```OFF_na``` Observed number of nonsilent mutations OFF target

- ```OFF_NA``` Number of nonsilent sites (corrected) OFF target

- ```OFF_ns``` Number of silent sites (corrected) OFF target

- ```OFF_NS``` Number of silent sites (corrected) OFF target


#### Notes before running
- No header needed for input VEP file
- VEP annotated first column must be in the format (chr_pos_ref/alt as described in https://www.ensembl.org/info/docs/tools/vep/vep_formats.html)
- VEP annotated file must only have chromosomes that are 1,2,3,4...22 or uppercase X,Y
- After you run vep with the option for ExAC frequencies, it would be necessary to remove all variants present in more than 0.1 percent of the population. You could apply the   filter using:
filter_vep -i input.annotated -f "ExAC_AF < 0.1 or not ExAC_AF" --ontology --filter "Consequence is coding_sequence_variant" 
-Be sure that you are using the GNU command line if you are running in a MacOS 
- Add dependencies to your path for easy running or hardcode the scripts
- The genome file used in SOPRANO or SSB is a two column file that contains the info of the name of the fasta id (column 1) and the length of that sequence (column 2).
- The UNIX system used should be able to recognize \t as a tab separator, some encodings may have problems on recognizing special characters
- Earlier versions of bedtools will not work
- Tab encoding should be \t (might be a problem for windows/OSX versions)
- Genome length file is a two column file specifying the fasta id and the length of the sequence (see how to obtain it at the bottom)
- Restrict your input dataset to chromosomes 1-22 and X and Y. Remove the rest.
- Input chromosome number must coincide with reference genome and annotation (chr1 vs 1)

[//]: # (## Genomes &#40;Updated&#41;)

[//]: # (To get hg19 &#40;GRCh37&#41; fasta genome, you can download it from ensembl. Now we provide the matching chrom sizes file in the data directory:)

[//]: # ()
[//]: # (```{bash})

[//]: # (wget https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.toplevel.fa.gz)

[//]: # ()
[//]: # (GENOME=~/tools/SOPRANO/data/chrom_GRCh37.sizes)

[//]: # (```)

## Obtain patient specific immune dN/dS values
To determine the patient specific immunopeptidome you should run the script 
`src/SOPRANO/scripts/get_epitope_HLA.pl`:

For example,
```{bash}
perl scripts/get_epitope_HLA.pl examples/TCGA_hlaTypesTEST.tsv 
```

This will create a command that you can run in a HPC cluster.
In our example case, the command is:

```{bash}
egrep -w -e "HLA-A2402|HLA-A0301|HLA-B1501|HLA-B1801|HLA-C0701|HLA-C0303|" data/allhlaBinders_exprmean1.IEDBpeps.mgd.bed | sortBed -i stdin | mergeBed -i stdin > TCGA-HQ-A5ND.exprmean1.IEDBpeps.SB.epitope.bed
egrep -w -e "HLA-A0205|HLA-A3303|HLA-B5301|HLA-B5301|HLA-C0401|HLA-C0401|" data/allhlaBinders_exprmean1.IEDBpeps.mgd.bed | sortBed -i stdin | mergeBed -i stdin > TCGA-FD-A6TC.exprmean1.IEDBpeps.SB.epitope.bed
```

In this previous case, the global immunopeptidome (containing all possible HLA 
alleles) was filtered by expression and similarity to positive T-cell assays 
epitopes from IEDB as described in the manuscript. Alternatively, users can 
obtain an unfiltered global immunopeptidome 
(now provided in https://www.synapse.org/#!Synapse:syn30371735) and perform 
the post-processing of each patient specific immunopeptidome with user-defined 
filters.  

After obtaining the immunopeptidome file following the steps before, you can 
run SOPRANO using the command following:

## Example SOPRANO pipeline run
```{bash}
cd src/SOPRANO
mkdir $HOME/test_run
RUN_SOPRANO -i examples/TCGA-05-4396-01A-21D-1855-08.annotated -b immunopeptidomes/human/TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed -n TCGA-05-4396 -o $HOME/test_run --use_ssb192 --exclude_drivers
```

## Limitations
- SOPRANO will fail if there are 0 synonymous mutations inside the immunopeptidome.
- Results must be taken cautiously when low number of synonymous mutations inside the immunopeptidome. See tutorial where we show the dependency of synonymous events versus the confidence interval length. In brief, one synonymous mutation in the immunopeptidome leads to a confidence interval length of ~0.6 and an 85% of misclassified cases
- Only mutations in those transcript IDs present on the immunopeptidome file will be used to estimate ON and OFF dN/dS. The rest of mutations will be discarded.
- SOPRANO can use mutations in intronic regions when using whole genome sequencing data to improve the estimate of the neutral mutation rate of the gene. However, this strategy is experimental since not every study will have information on all possible mutations falling into intronic regions (the length of the intronic region was estimated using hg19 reference). 

