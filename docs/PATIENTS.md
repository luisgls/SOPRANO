## Obtain patient specific immune dN/dS values
To determine the patient specific immunopeptidome you should run the script 
`src/SOPRANO/scripts/get_epitope_HLA.pl`:

For example,
```{bash}
perl scripts/get_epitope_HLA.pl examples/TCGA_hlaTypesTEST.tsv 
```

This will create a command that you can run in a HPC cluster.
In our example case, the command is:

```{bash}
egrep -w -e "HLA-A2402|HLA-A0301|HLA-B1501|HLA-B1801|HLA-C0701|HLA-C0303|" data/allhlaBinders_exprmean1.IEDBpeps.mgd.bed | sortBed -i stdin | mergeBed -i stdin > TCGA-HQ-A5ND.exprmean1.IEDBpeps.SB.epitope.bed
egrep -w -e "HLA-A0205|HLA-A3303|HLA-B5301|HLA-B5301|HLA-C0401|HLA-C0401|" data/allhlaBinders_exprmean1.IEDBpeps.mgd.bed | sortBed -i stdin | mergeBed -i stdin > TCGA-FD-A6TC.exprmean1.IEDBpeps.SB.epitope.bed
```

In this previous case, the global immunopeptidome (containing all possible HLA 
alleles) was filtered by expression and similarity to positive T-cell assays 
epitopes from IEDB as described in the manuscript. Alternatively, users can 
obtain an unfiltered global immunopeptidome 
(now provided in https://www.synapse.org/#!Synapse:syn30371735) and perform 
the post-processing of each patient specific immunopeptidome with user-defined 
filters.  

After obtaining the immunopeptidome file following the steps before, you can 
run SOPRANO using the command following:

## Example SOPRANO pipeline run
```{bash}
cd src/SOPRANO
mkdir $HOME/test_run
RUN_SOPRANO -i examples/TCGA-05-4396-01A-21D-1855-08.annotated -b immunopeptidomes/human/TCGA-05-4396.Expressed.IEDBpeps.SB.epitope.bed -n TCGA-05-4396 -o $HOME/test_run --use_ssb192 --exclude_drivers
```