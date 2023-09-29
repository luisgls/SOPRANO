# Downloading the genome reference files

SOPRANO has a built-in genome reference downloader for Homo Sapiens 
(GRCh37 and GRCh38). To download a reference genome, run the command:

```GET_GENOMES -r <reference>```

Currently supported references are `GRCh37` and `GRCh38`. By default, the 
download will pull the ensembl 110 data release. Users can prescribe a 
different release if they wish via the `--release` flag:

```GET_GENOMES -r GRCh37 --release 109```

Users will be prompted to download the toplevel and primary assembly genome
files. 


### Notes


- **Toplevel genome references are a core dependency of SOPRANO pipelines**, 
whilst primary assemblies are simply used to accelerate the annotation 
procedure (currently in development as of September 2023).


- If users would like to link their existing VEP cache containing reference
genome files, see the accompanying [documentation](VEP.md).


- Once references are downloaded, SOPRANO will automatically compute the 
corresponding fasta index files and chromosome sizes. Note that since these 
files are lightweight, some are readily shipped with the SOPRANO pacakge.


- The downloaded data mimics the cache structure of ensembl VEP, i.e., files will
be downloaded to `/src/SOPRANO/data/homo_sapiens/<release>_<reference>/`.


- **Warning:** Downloading all the files associated with a given reference will
consume ~70GB of memory. Users are therefore advised to only download the 
primary assembly files if they are required for annotations.