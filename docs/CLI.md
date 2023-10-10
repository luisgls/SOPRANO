# Command line interface

SOPRANO can be run from the command line, following the syntax:

```{shell}
soprano-run <arguments>
```

### Required arguments

These arguments have short (`-s`) or long (`--long`) argument flags.

- `-n | --name`

    The name of the SOPRANO pipeline run.


- `-o | --output`

    The path to the output file directory (i.e. the cache) for the SOPRANO
    pipeline run.


- `-i | --input`

    The path to the input VEP annotated file that defines mutations.


- `-b | --bed_file`

    The path to the bed file containing protein coordinates named by their 
    transcript ID (e.g. ENSTXXXXXX 123 135)


### Optional arguments

These arguments have short (`-s`) or long (`--long`) argument flags.

- `-m | --random_regions`

    The path to the bed file that defines regions for randomization.
    Defaults to None, implying arbitrary randomization.


- `-s | --seed`

    The seed value to be used in randomization process. If defined, users 
    should be able to reproduce their results when randomization is applied.
    This value should be an integer greater than zero. Defaults to "-1", 
    implying no seed value chosen.

  
- `-s | --species`

    The Latin species name associated with genome reference. The default value
    is for homo_sapiens.


- `-a | --assembly`

    The name of the Ensembl reference genome assembly. The default value is for 
    GRCh38 (Homo Sapiens).


- `-r | --release`

    The Ensembl release number of the genome reference
    file. The default value is 110.


- `-t | --transcript`

    Path to transcript lengths.
    ([Default](../src/SOPRANO/data/transcript_length.txt))


- `-p | --protein_transcript`

    Path to protein transcript lengths.
    ([Default](../src/SOPRANO/data/protein_length.txt))


- `-f | --fasta`

  Path to ensembl transcript IDs (fasta file).
  ([Default](../src/SOPRANO/data/ensemble_transcriptID.fasta))


### Optional switches

These switches are activated by passing the long (`--long`) flag to the command
line with no subsequent argument. They alter the runtime behaviour of SOPRANO.

- `--use_random`

    Run SOPRANO with randomization. Note that if `random_regions` are passed,
    to the CLI, this overrides the selection to `use_random` to `True`.


- `--use_ssb192`

    Run SOPRANO using SSB192 substitutions, otherwise defaults SSB7.


- `--keep_drivers`

    Run SOPRANO with driver genes included. Otherwise, 
    [driver genes](../src/SOPRANO/data/genes2exclude.txt) are excluded.