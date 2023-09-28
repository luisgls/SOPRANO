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