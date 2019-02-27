#!/usr/bin/perl

open(FH, $ARGV[0]);
use strict;

#my $DIR="~/Projects/Immunopeptidome/binder";
my $database="data/allhlaBinders_exprmean1.IEDBpeps.mgd.bed";
#my $database="/home/haem/lortiz/Projects/Immunopeptidome/run_all/allhlaBinders.IEDBpeps.bed";
my $query="";

while (my $line = <FH>) {
    chomp $line;
    my ($id,$ttype,$hla)= split /\t/ , $line;
    my @hlas=split(',', $hla);
    #print Dumper(@hlas);
    foreach my $i (@hlas){
    	chomp $i;
	$i =~ s/\r//g;
	#print $i."\n";
	$query .= "$i|";
	}
	print "egrep -w -e \"";
	print $query;

    print "\" $database | sortBed -i stdin | mergeBed -i stdin > $id.exprmean1.IEDBpeps.SB.epitope.bed\n";  
    #print("cat $DIR/$i.noMGD.bed.expmean1med1.bed  | sortBed -i stdin | mergeBed -i stdin >> $id.epitope.expfilt.bed\n");   
    $query="";	
    #print("sortBed -i $id.epitope.bed | mergeBed -i stdin > $id.mgd.bed\n");		
#	print("sortBed -i $id.epitope.expfilt.bed | mergeBed -i stdin > $id.mgd.expfilt.bed\n");		
}
