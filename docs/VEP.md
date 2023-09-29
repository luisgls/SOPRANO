# Linking existing VEP cache(s)

SOPRANO is configured to internally search for genome reference files inside
`./src/SOPRANO/data`. If SOPRANO is used to [download references](GENOMES.md),
then this is where they will be downloaded to, following the standard VEP
cache pattern `*/<species>/<ref_release>/<files>`.

If users of SOPRANO are existing users of VEP, it is possible to link the fasta
and fasta index files to the SOPRANO `data` folder.

From the command line, simply run the command

```shell
soprano-link-vep
```

which will link files found inside of `$HOME/.vep` and compute the chromosome
sizes from the toplevel reference file (if found). If users have a VEP cache
other than the default `$HOME/.vep` location, this can be linked via specifying

```shell
soprano-link-vep --cache /path/to/vep/cache
```