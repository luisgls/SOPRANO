#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my %freq1 = GetFreq1($ARGV[0]) or die "missing per tumor frequencies";

my %freq2 = GetFreq2($ARGV[1]) or die "missing freq per mut type in each ttype";

#print Dumper(%freq2);

my $file=$ARGV[0];

open(FH, $file) or die "Not able to open file";

while (my $line = <FH>) {
    chomp $line;
    my ($info,$oriNa,$oriNs)=split("\t",$line);
    my (undef,$change,$compl)=split("_",$info);
    ##For each tumor type multiplication-context factor
    my $rateparam=$freq2{$change}{val}/$freq1{$change}{val};
    
    print join("\t",$change,$rateparam,$oriNa,$oriNa*$rateparam,$oriNs,$oriNs*$rateparam)."\n";
    
}


sub GetFreq1{
    my $file=shift;
    my %freq1;
    my $total=0;
    open( FQ, $file );
    while (my $line = <FQ>){
        next if $line =~ m/\#/;
        chomp $line;
        my @tmp = split("\t|_",$line);
       ## get frequency of mutations per tumor
        $freq1{$tmp[1]}{ind}=$tmp[3]+$tmp[4];
        $total += $freq1{$tmp[1]}{ind};
    }
    
    foreach my $final(keys %freq1){
        $freq1{$final}{val}=$freq1{$final}{ind}/$total;
    }
    return %freq1;
   
}

sub GetFreq2{
    my $file=shift;
    my %freq2;
    my $total=0;
    open( FQ, $file );
    while (my $line = <FQ>){
        next if $line =~ m/\#/;
        chomp $line;
        my @tmp = split("\t",$line);
       ## get frequency of mutations per tumor
        $freq2{$tmp[1]}{ind}=$tmp[2];
        $total += $freq2{$tmp[1]}{ind};
    }
    
    foreach my $final(keys %freq2){
        $freq2{$final}{val}=$freq2{$final}{ind}/$total;
    }
    return %freq2;
}




sub getNs{
    my $file=shift;
    my %info;
    open(FS, $file);
    while (my $line2 = <FS>) {
        chomp $line2;
        my @tmp2=split("\t",$line2);
        
        $info{$tmp2[0]}{Na}=$tmp2[1];
        $info{$tmp2[0]}{Ns}=$tmp2[2];
    }
    return %info;
}


sub getTypeFreq{
    my $file=shift;
    my %ttype;
    open(FA, $file);
    while (my $line3 = <FA>) {
        chomp $line3;
        my ($tumor,$chng,$fold)=split("\t",$line3);
        $ttype{$tumor}{$chng}=$fold;
    }
    return %ttype;
}
