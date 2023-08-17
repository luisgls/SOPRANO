#! /usr/bin/perl
use strict;
use Data::Dumper;

my %hash_ttype; 
#my (%DNA_complement) = ( 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'cpg' => 'rev');
my (%variations) = ( 
			'A' =>'TCG',
			'T' => 'ACG',
			'C' => 'TAG',
			'G' => 'ATC' 
			);
my (%genetic_code) = (
			'TCA' => 'S', # Serine Polar
			'TCC' => 'S', # Serine Polar
			'TCG' => 'S', # Serine Polar
			'TCT' => 'S', # Serine Polar
			'TTC' => 'F', # Phenylalanine HP
			'TTT' => 'F', # Phenylalanine HP
			'TTA' => 'L', # LeTcine HP
			'TTG' => 'L', # LeTcine HP
			'TAC' => 'Y', # Tyrosine Polar
			'TAT' => 'Y', # Tyrosine Polar
			'TAA' => '*', # Stop
			'TAG' => '*', # Stop
			'TGC' => 'C', # Cysteine Polar
			'TGT' => 'C', # Cysteine Polar
			'TGA' => '*', # Stop
			'TGG' => 'W', # Tryptophan HP
			'CTA' => 'L', # LeTcine HP
			'CTC' => 'L', # LeTcine HP
			'CTG' => 'L', # LeTcine HP
			'CTT' => 'L', # LeTcine HP
			'CCA' => 'P', # Proline HP
			'CCC' => 'P', # Proline HP
			'CCG' => 'P', # Proline HP
			'CCT' => 'P', # Proline HP
			'CAT' => 'H', # Histidine Charged
			'CAC' => 'H', # Histidine Charged
			'CAA' => 'Q', # GlTtamine Charged
			'CAG' => 'Q', # GlTtamine Charged
			'CGA' => 'R', # Arginine Charged
			'CGC' => 'R', # Arginine Charged
			'CGG' => 'R', # Arginine Charged
			'CGT' => 'R', # Arginine Charged
			'ATA' => 'I', # IsoleTcine HP
			'ATC' => 'I', # IsoleTcine HP
			'ATT' => 'I', # IsoleTcine HP
			'ATG' => 'M', # Methionine HP
			'ACA' => 'T', # Threonine Polar
			'ACC' => 'T', # Threonine Polar
			'ACG' => 'T', # Threonine Polar
			'ACT' => 'T', # Threonine Polar
			'AAC' => 'N', # Asparagine Polar
			'AAT' => 'N', # Asparagine Polar
			'AAA' => 'K', # Lysine Charged
			'AAG' => 'K', # Lysine Charged
			'AGC' => 'S', # Serine Polar
			'AGT' => 'S', # Serine Polar
			'AGA' => 'R', # Arginine Charged
			'AGG' => 'R', # Arginine Charged
			'GTA' => 'V', # Valine HP
			'GTC' => 'V', # Valine HP
			'GTG' => 'V', # Valine HP
			'GTT' => 'V', # Valine HP
			'GCA' => 'A', # Alanine HP
			'GCC' => 'A', # Alanine HP
			'GCG' => 'A', # Alanine HP
			'GCT' => 'A', # Alanine HP
			'GAC' => 'D', # Aspartic Acid 
			'GAT' => 'D', # Aspartic Acid
			'GAA' => 'E', # GlTtamic Acid
			'GAG' => 'E', # GlTtamic Acid
			'GGA' => 'G', # Glycine HP
			'GGC' => 'G', # Glycine HP
			'GGG' => 'G', # Glycine HP
			'GGT' => 'G'  # Glycine HP
			);
	
	##### WE WILL GET ALL INPUT FASTA SEQUENCES IN FRAME#####
	open (HUMAN, $ARGV[0]) or die; 
	my %sequences;
	while (<HUMAN>)
	{
		my $read = $_;
		chomp $read;
		if ($read =~ m/>/)
		{
			
			my ($id) = ($read =~ m/>(\S+)/);
			$read = <HUMAN>;
			chomp $read;
			my $seq = $read;
			$sequences{$id} = $seq;
		}
		
	}
	

	
	###### First step --> loop through all transcripts found in that type of cancer ###########
	open (TRANSCRIPT_LIST, $ARGV[1]) or die "No transcript or fasta file provided";
	while (<TRANSCRIPT_LIST>)
	{
		my $read = $_; 
		chomp $read;
		
		my @sequence_trans = split ('', $sequences{$read});
		my @normal_positions;
		for my $i (0 .. $#sequence_trans)
		{
	      		push (@normal_positions, $i+1);
	      	}
		#print Dumper(@sequence_trans);
		
		my @aa = ($sequences{$read}=~ m/.../g );
		#print Dumper(@aa);
		
		###### LETS GET ALL THE SITES FOR EACH SPECIFIC VARIATION ####
		my %sites_ratio_syn;
		my %sites_ratio_nonsyn;
		my %tricontext;
		my $syn;
		my $nonsyn;
		
		#my $sequence_pos = 0;
		my $possible1;
		my $possible2;
		my $possible3;
		my @all_possible1;
		my @all_possible2;
		my @all_possible3;
		my $triplet_count = 0;
		my $nucleotide_count = 0;
		#print Dumper($sequence_trans[$nucleotide_count]);
		foreach my $triplet (@aa)
		{

			my @list = split ('', $triplet);
			my $nucl1 = $list[0]; #first element of the triplet
			my $nucl2 = $list[1]; #second element of the triplet
			my $nucl3 = $list[2]; #third element of the triplet
			
			$nucleotide_count++;

			$possible1 = $variations{$nucl1}; #all possible mutations for the first nucleotide(so A to C,G,T and so on)
			@all_possible1 = split ('', $possible1);
			foreach my $thing (@all_possible1)
			{
				my $actual_change = join ('/', $nucl1, $thing);
				my $new_codon = join ('', $thing, $nucl2, $nucl3);
				my $final_change = join("_",$triplet,$new_codon);

				next if ($nucleotide_count == 1);

				my $context=join("",$sequence_trans[$nucleotide_count-2],$sequence_trans[$nucleotide_count-1],$sequence_trans[$nucleotide_count]."_".$thing);
				$tricontext{$context}=1;
				next if $triplet_count == 0;
				next if $triplet_count == scalar(@aa)-1;
					
				if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
				{
					$sites_ratio_nonsyn{$context} = $sites_ratio_nonsyn{$context} +1;
				}
				else
				{
					$sites_ratio_syn{$context} = $sites_ratio_syn{$context} +1;
				}
			}

			$nucleotide_count++;

			$possible2 = $variations{$nucl2}; #all possible mutations for the second nucleotide(so A to C,G,T and so on)
			@all_possible2 = split ('', $possible2);
			foreach my $thing (@all_possible2)
			{
				my $actual_change = join ('/', $nucl2, $thing);
				my $new_codon = join ('', $nucl1, $thing, $nucl3);
				my $final_change = join("_",$triplet,$new_codon);
				my $context=join("",$sequence_trans[$nucleotide_count-2],$sequence_trans[$nucleotide_count-1],$sequence_trans[$nucleotide_count]."_".$thing);
				$tricontext{$context}=1;
				next if $triplet_count == 0;
				next if $triplet_count == scalar(@aa)-1;
				#print Dumper($triplet_count);
				#print Dumper(scalar(@aa));
				if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
				{
					$sites_ratio_nonsyn{$context}++;
				}
				else
				{
					$sites_ratio_syn{$context}++;
				}
			}

		$nucleotide_count++;

		$possible3 = $variations{$nucl3};  #all possible mutations for the third nucleotide(so A to C,G,T and so on)
		@all_possible3 = split ('', $possible3);
			foreach my $thing (@all_possible3)
				{
				my $actual_change = join ('/', $nucl3, $thing);
				my $new_codon = join ('',  $nucl1, $nucl2, $thing);
				my $final_change = join("_",$triplet,$new_codon);
				next if ($nucleotide_count == scalar(@sequence_trans));
				my $context=join("",$sequence_trans[$nucleotide_count-2],$sequence_trans[$nucleotide_count-1],$sequence_trans[$nucleotide_count]."_".$thing);
				$tricontext{$context}=1;
				next if $triplet_count == 0;
				next if $triplet_count == scalar(@aa)-1;

				if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
				{
					$sites_ratio_nonsyn{$context}++;
				}
				else
				{
					$sites_ratio_syn{$context}++;
				}
			}
			$triplet_count++;
		}
		my %team_x = (%sites_ratio_nonsyn, %sites_ratio_syn);	
		my $nonsyn_sites_total=0;
		my $syn_sites_total=0;
		foreach my $key(sort keys %team_x)
		{
			my $nonsyn_sites_total = ($sites_ratio_nonsyn{$key})/3;
			my $syn_sites_total = ($sites_ratio_syn{$key})/3;
			if ($key =~ m/A\/C|A\/G|A\/T|C\/A|C\/G|C\/T|cpg/){
				print "$read\t$key\t$key\t$nonsyn_sites_total\t$syn_sites_total\n"; 
			}else{
				print "$read\t$key\t$key\t$nonsyn_sites_total\t$syn_sites_total\n"; 
			}				
		}
	}

