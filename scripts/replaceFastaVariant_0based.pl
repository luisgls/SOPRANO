#!/usr/bin/perl
use Bio::Seq;
use Bio::SeqIO;
use Data::Dumper;
use strict;

open(FH, $ARGV[0]);

my %seq=getSeq($ARGV[1]);

my %countMuts=countM($ARGV[0]);

my %bed=getTarget($ARGV[2]);
#print Dumper(%bed);


#Read orinially annotated VEP file
while (my $line = <FH>){
    chomp $line;
    my ($vepid,$vepid2,$alt,$geneid,$enstid,undef,$vartype,$posENSG,$posENST,$posENSP,$aarefalt,$codonrefalt,undef,$info)=split("\t",$line);   
    my ($refaa,$altaa)= split("/",$aarefalt);
    
    my $protseqstart="";
    my $protseqend="";
    my $finalseq="";
    
    next if ($line =~ m/^#/);
    
    if ($posENSP > 0 && length($aarefalt) > 1){
        $protseqstart   = substr($seq{$enstid}{'seq'},0,$posENSP-1).$altaa;
        $protseqend     = substr($seq{$enstid}{'seq'},$posENSP,length($seq{$enstid}{'seq'}));
        $finalseq       = $protseqstart.$protseqend;
      }
    else {
        #print $line."\t"."REFERENCE"."\n";
        next;
    }
    
    

    if (exists $bed{$enstid}{in}) {
        foreach my $dump( @{$bed{$enstid}{coord}}){
            my ($start,$stop)=split("_",$dump);    
            my ($targetlengthwt)=($stop-$start)+1;
            my ($targetlengthmut)=($stop-$start)+1;
            
            if ($altaa eq "*") {
                $targetlengthmut=($posENSP-$start)+1;
            }
            
            my $wtPEP=substr($seq{$enstid}{'seq'},$start-1,$targetlengthwt);
            my $mutPEP=substr($finalseq,$start-1,$targetlengthmut);
            
            
            if ($posENSP>=$start && $posENSP <= $stop) {
                #print $line."\t".$finalseq."\t".substr($seq{$enstid}{'seq'},$start-1,$targetlength)."\t".substr($finalseq,$start-1,$targetlength)."\t".$start."-".$stop."\n";
                print $line."\t".$wtPEP."\t".$mutPEP."\t".$start."-".$stop."\n";
                #last;
            }else{
                next;   
            }           
        }
    }
    else{
         #print $line."\t".$finalseq."\n";
    }        
}



sub getSeq{
my $in = Bio::SeqIO->new('-file' => $_[0], '-format' => 'Fasta');
my %hash;
while (my $seq = $in->next_seq()) {
     $hash{$seq->display_id()}{'seq'}=$seq->seq;
}
return %hash;
}


sub countM{
    open(FH2,shift);
    my %hash;
    while (my $line2 = <FH2>) {
        chomp $line2;
        my @tmp=split("\t",$line2);   
        $hash{$tmp[4]}{count}++;
    }
    return %hash;
    
}

sub getTarget{
    open(FH3,shift);
    my %hashis;
    while (my $line3 = <FH3>) {
        chomp $line3;
        my @tmp=split("\t",$line3);
        $hashis{$tmp[0]}{in}=1;
        
        push (@{$hashis{$tmp[0]}{coord}},$tmp[1]."_".$tmp[2]);
    }
    return %hashis;
}
