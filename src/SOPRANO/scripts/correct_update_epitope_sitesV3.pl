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
    my (undef,$change1,$compl,undef,undef)=split("_",$info);
    my $change=join("_",$change1,$compl);
    ##For each tumor type multiplication-context factor
    #print Dumper(%freq2);
    #print Dumper($freq2{$change}{val});
    my $rateparam=1;
    if (exists $freq2{$change}{val}){
        $rateparam=($freq2{$change}{val}+1/$freq2{total})/($freq1{$change}{val}+1/$freq1{total});
    }
    else{
        # When trinuclecotide file doesnt have counts of ttrinuc changes assume lowest frequency.
        # print STDERR "Frequency of $change is 0, setting to the lowest frequency 1/$freq2{total}\n ";
        $rateparam=(1/$freq2{total})/($freq1{$change}{val}+1/$freq1{total});
    }
    
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
        $freq1{$tmp[1]."_".$tmp[2]}{ind}=$tmp[5]+$tmp[6];
        $total += $freq1{$tmp[1]."_".$tmp[2]}{ind};
    }
    
    foreach my $final(keys %freq1){
        $freq1{$final}{val}=$freq1{$final}{ind}/$total;
        $freq1{total}=$total;
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
        $freq2{$tmp[0]}{ind}=$tmp[1];
        $total += $freq2{$tmp[0]}{ind};
    }
    
    foreach my $final(keys %freq2){
        $freq2{$final}{val}=$freq2{$final}{ind}/$total;
        $freq2{total}=$total;
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
