#!/usr/bin/perl
use strict;

use warnings;

my $file=$ARGV[0];

my %hash = getInfo($ARGV[1]);

my @tmp;

#print Dumper(%hash);

open DOT, $file;
while (my $line = <DOT>) {
    
    chomp $line;
    my @tmp = split("\t",$line);
    
    my $context=$tmp[0];
    #my ($ref,$alt) = split("/",$tmp[1]);
    my $base = $context;
    
    ##print $tmp[1];
    if(exists $hash{$base}){
        print join("\t",@tmp,$hash{$base});
        }
    else{
        print "WARNING: base not found\n";
        print join("\t",@tmp,$tmp[1])."\n";
    }
}

sub getInfo{
    my %hash;
    open FH , shift;
    
    while (my $info = <FH>){
            my ($map,$map2)=split("\t",$info);
            $hash{$map}=$map2;
    }
    return %hash;        
}


