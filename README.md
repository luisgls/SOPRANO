# SOPRANO: Selection On PRotein ANnotated regiOns
The SOPRANO method was developed to quantify selection in specific regions of 
the genome (Zapata et al, in revision). It calculates ON- and OFF-target dN/dS 
using a set of annotated somatic point mutations and a transcript coordinates 
file. 

_What inputs does SOPRANO require?_

- A set of mutations (missense/truncating) and their respective functional 
annotation (derived from ensembl VEP). Mutations can be filtered a priori by the user 
(i.e. by predicted binding affinity or by expression status).

- A set of ensembl transcript coordinates where selection will be estimated.

_How does SOPRANO calculate dN/dS?_

- "ON" dN/dS is the value calculated inside the coordinates provided using a 
192-trinucleotide correction signature obtained "on-the-fly" from the input 
mutation file. Alternatively, the user can provide a pre-calculated 
trinucleotide mutation frequency file. Importantly, ON dN/dS and OFF dN/dS 
(the portion outside the coordinates provided) will be calculated only using 
transcripts defined in this file. 

SOPRANO can be run from the command line, or via a streamlit app interface. 
There are additional tools built within the SOPRANO installation that enable 
users to download genomes, link existing VEP caches, and annotate VCF files.

## Documentation

- [Installation](docs/INSTALL.md) - Building the environment and installing SOPRANO
- [CLI](docs/CLI.md) - Running SOPRANO from the command line.
- [GUI](docs/APP.md) - Running SOPRANO from your web browser.
- [Genomes](docs/GENOMES.md) - Downloading and preparing reference genome files.
- [Link VEP](docs/VEP.md) - Linking your existing VEP cache files.
- [Inputs 1](docs/INPUT.md) - Overview of the VEP annotated input.
- [Inputs 2](docs/BED.md) - Overview of the immunopeptidome input.
- [Output](docs/OUTPUT.md) - Overview of the summary statistics computed by SOPRANO.
- [Patient Specific](docs/PATIENTS.md) - How to obtain a patient specific immune dN/dS.
- [Notes](docs/NOTES.md) - Notes on and limitations of this software.

### Contact

- [Luis Zapata](mailto:Luis.Zapata@icr.ac.uk)
- [Kareem Marzouk](mailto:Kareem.Marzouk@icr.ac.uk) / [ICR Scientific Software Group](mailto:scientificcomputingteam@icr.ac.uk)