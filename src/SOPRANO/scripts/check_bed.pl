#!/usr/bin/perl -w
use strict;

open(FH,$ARGV[0]) or die "input file missing";

while (my $line = <FH>) {
    chomp $line;
    my @tmp = split ("\t",$line);
    
    if (scalar @tmp < 4) {
        print $line."\nERROR:Wrong bed number of columns\n";
    }
    
    
    #Check for chromosome
    #if ($tmp[0] =~ m/[1-22]|[XY]/) {
    #    print "Chromosome OK\n";
    #}else{
    #    print $line."\nChromosome integer wrong in vep.annotated file\n";
    #}
    
    #check for coordinates start 
=cut
    if($tmp[1]>$tmp[2]){
	print join("\t",$tmp[0],$tmp[2],$tmp[1]);
        splice(@tmp,0,3);
	foreach my $l(@tmp){
	    print "\t".$l;
		}
	print "\n";
	}
	else{
	    foreach my $l(@tmp){
	    print $l."\t";
		}
	print "\n";
	}
=cut
    
    
    #check for wrong characters
    
    if ($tmp[1] =~ m/\D/i || $tmp[2] =~ m/\D/i) {
        print $line."\nERROR2:Bed coordinates not integers Check\n";
    }else{
        #print "Bed coordinates only integers OK\n";
    }
    
}
    