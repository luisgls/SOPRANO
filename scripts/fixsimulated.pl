#!/usr/bin/perl -w
use strict;


open(FH, $ARGV[0]);


while (my $line = <FH>) {
    chomp $line;
    
    if ($line =~ m/#/){
        print $line."\n";
        next;
    }
    
    
    my @tmp = split("\t",$line);
    
    if ($tmp[13] =~ m/STRAND=-1/) {
        $tmp[0] =~ tr/TCGA/AGCT/;
        $tmp[2] =~ tr/TCGA/AGCT/;
    }
    
    
    my (@info)=split(":",$tmp[0]);
    my ($chr,$pos)=split(":",$tmp[1]);
    
    
    if (scalar @info == 2) {
        my $change = $info[1];
        $change =~ s/c.[0-9]+//g;
        $change =~ s/>/\//;
        shift @tmp for 1..1;
        print join("\t",$chr."_".$pos."_".$change,@tmp)."\n";
    }
    else{
        print $line."\n";
    }
    
}

