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