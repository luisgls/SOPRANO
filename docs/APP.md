# Graphical application 

SOPRANO can be run from your local browser via a streamlit app, simply run

```{shell}
soprano-app
```


### Running the SOPRANO pipeline

The App interface searches for genome reference files, annotated input files,
and BED protein files readily available on the users system, and presents
them as drop down selections to the user. 

It is straightforward to specify the run time configuration of SOPRANO 
(e.g., to exclude driver genes) with tick boxes. Once a configuration is 
selected, users can execute the SOPRANO pipeline via the **Run Pipeline** 
button.

The processes deployed via the app are identical to those implemented
via the [CLI](CLI.md) - the pipeline steps are equivalent.


### Pipeline data cache

By default, the application will cache pipeline data into the folder relative 
to the SOPRANO repository root named `./pipeline_cache`.

This selection can be overriden by setting the environment variable 
`SOPRANO_CACHE` prior to launching the app. For example,

```shell
export SOPRANO_CACHE=/path/to/my/soprano/cache
soprano-app
```

### Data sources

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

This app is actively being developed (September 2023). The configuration is 
subject to change, and additional features are soon to be integrated.

_Planned features:_

- Generalization of species. Currently only Homo Sapiens are supported.
- Selection of transcript files. Currently only available via CLI.
- Interface for linking VEP cache files.
- Interface for downloading genome references.

### Example app

![SOPRANO App](soprano_app.PNG "SOPRANO App")