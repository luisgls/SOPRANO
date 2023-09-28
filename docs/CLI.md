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

    The name of the genome reference. Currently, this CLI supports GRCh37 and 
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