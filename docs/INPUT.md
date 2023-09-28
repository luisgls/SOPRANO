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