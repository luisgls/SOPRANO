#### Notes before running
- No header needed for input VEP file
- VEP annotated first column must be in the format (chr_pos_ref/alt as described in https://www.ensembl.org/info/docs/tools/vep/vep_formats.html)
- VEP annotated file must only have chromosomes that are 1,2,3,4...22 or uppercase X,Y
- After you run vep with the option for ExAC frequencies, it would be necessary to remove all variants present in more than 0.1 percent of the population. You could apply the   filter using:
filter_vep -i input.annotated -f "ExAC_AF < 0.1 or not ExAC_AF" --ontology --filter "Consequence is coding_sequence_variant" 
-Be sure that you are using the GNU command line if you are running in a MacOS 
- Add dependencies to your path for easy running or hardcode the scripts
- The genome file used in SOPRANO or SSB is a two column file that contains the info of the name of the fasta id (column 1) and the length of that sequence (column 2).
- The UNIX system used should be able to recognize \t as a tab separator, some encodings may have problems on recognizing special characters
- Earlier versions of bedtools will not work
- Tab encoding should be \t (might be a problem for windows/OSX versions)
- Genome length file is a two column file specifying the fasta id and the length of the sequence (see how to obtain it at the bottom)
- Restrict your input dataset to chromosomes 1-22 and X and Y. Remove the rest.
- Input chromosome number must coincide with reference genome and annotation (chr1 vs 1)



## Limitations
- SOPRANO will fail if there are 0 synonymous mutations inside the immunopeptidome.
- Results must be taken cautiously when low number of synonymous mutations inside the immunopeptidome. See tutorial where we show the dependency of synonymous events versus the confidence interval length. In brief, one synonymous mutation in the immunopeptidome leads to a confidence interval length of ~0.6 and an 85% of misclassified cases
- Only mutations in those transcript IDs present on the immunopeptidome file will be used to estimate ON and OFF dN/dS. The rest of mutations will be discarded.
- SOPRANO can use mutations in intronic regions when using whole genome sequencing data to improve the estimate of the neutral mutation rate of the gene. However, this strategy is experimental since not every study will have information on all possible mutations falling into intronic regions (the length of the intronic region was estimated using hg19 reference). 

