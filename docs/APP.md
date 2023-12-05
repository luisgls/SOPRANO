# Graphical application

SOPRANO can be run from your local browser via a streamlit app, simply run

```shell
soprano-app
```

The SOPRANO interface is organized into a collection of steps

- Step 1: Download genome reference files
- Step 2: Annotate mutation sources using VCF input files
- Step 3: Prepare immunopeptidome file
- Step 4: Run pipeline

### App data sources

By default, the application will serve options for the annotated mutation
and human immunopeptidome files shipped with SOPRANO. Additional files
can be detected by the application by their placement in the `./app_sources`
folders:

- `./app_sources/annotated_inputs`

  VEP annotated mutation files placed in this directory will be detected,
  so long as they have the extension pattern `*anno*`. E.g., `mutations.anno`
  or `mutations.annotated`.


- `./app_sources/immunopeptidomes`

  User defined immunopeptidomes BED files placed in this directory will be
  detected, so long as they have the extension `.bed`. E.g., `immuno.bed`.


- `./app_sources/coordinate_files`

  User defined BED files that can be used for randomization will be detected,
  so long as they have the extension `.bed`. E.g., `randoms.bed`.

### Development status

This app is actively being developed (November 2023). The configuration is
subject to change, and additional features are soon to be integrated.