#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;


my %hash = GetFreq($ARGV[1]) or die "missing per tumor frequencies";
#my %bigN = getNs($ARGV[2]) or die "missing per gene big Ns values";
my %ttypefreq = getTypeFreq($ARGV[2]) or die "missing freq per mut type in each ttype";

#print Dumper(%ttypefreq);

my $file=$ARGV[0];

open(FH, $file) or die "Not able to open file";

while (my $line = <FH>) {
    chomp $line;
    my ($info,$oriNa,$oriNs)=split("\t",$line);
    my (undef,$change,$compl)=split("_",$info);
    ##For each tumor type multiplication-context factor
    foreach my $caca(sort keys %ttypefreq){
        
        if (exists $ttypefreq{$caca}{$change}){
            print $caca."\t".$change."\t".$ttypefreq{$caca}{$change}.
            "\t".$oriNa."\t".$oriNa*$ttypefreq{$caca}{$change}.
            "\t".$oriNs."\t".$oriNs*$ttypefreq{$caca}{$change}.
            "\t".$hash{$caca}{val}*$oriNa*$ttypefreq{$caca}{$change}.
            "\t".$hash{$caca}{val}*$oriNs*$ttypefreq{$caca}{$change}.
            "\n";
        }
        else{
            print "Values Not Found\n";
        }
        
    }
}
sub GetFreq{
    my $file=shift;
    my %freq;
    my $total=0;
    open( FQ, $file );
    while (my $line = <FQ>){
        next if $line =~ m/\#/;
        chomp $line;
        my @tmp = split("\t",$line);
       ## get frequency of mutations per tumor
        $freq{$tmp[0]}{ind}=$tmp[3];
        $total += $freq{$tmp[0]}{ind};
    }
    
    foreach my $final(keys %freq){
        $freq{$final}{val}=$freq{$final}{ind}/$total;
    }
    return %freq;
   
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
