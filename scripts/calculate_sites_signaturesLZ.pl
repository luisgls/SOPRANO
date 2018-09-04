#! /usr/bin/perl
use strict;
use Data::Dumper;

my %hash_ttype; 
my (%DNA_complement) = ( 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', 'cpg' => 'rev');
my (%variations) = ( 
			'A' =>'TCG',
			'T' => 'ACG',
			'C' => 'TAG',
			'G' => 'ATC' 
			);
my (%genetic_code) = (
			'TCA' => 'S', # Serine
			'TCC' => 'S', # Serine
			'TCG' => 'S', # Serine
			'TCT' => 'S', # Serine
			'TTC' => 'F', # Phenylalanine
			'TTT' => 'F', # Phenylalanine
			'TTA' => 'L', # LeTcine
			'TTG' => 'L', # LeTcine
			'TAC' => 'Y', # Tyrosine
			'TAT' => 'Y', # Tyrosine
			'TAA' => '*', # Stop
			'TAG' => '*', # Stop
			'TGC' => 'C', # Cysteine
			'TGT' => 'C', # Cysteine
			'TGA' => '*', # Stop
			'TGG' => 'W', # Tryptophan
			'CTA' => 'L', # LeTcine
			'CTC' => 'L', # LeTcine
			'CTG' => 'L', # LeTcine
			'CTT' => 'L', # LeTcine
			'CCA' => 'P', # Proline
			'CCC' => 'P', # Proline
			'CCG' => 'P', # Proline
			'CCT' => 'P', # Proline
			'CAT' => 'H', # Histidine
			'CAC' => 'H', # Histidine
			'CAA' => 'Q', # GlTtamine
			'CAG' => 'Q', # GlTtamine
			'CGA' => 'R', # Arginine
			'CGC' => 'R', # Arginine
			'CGG' => 'R', # Arginine
			'CGT' => 'R', # Arginine
			'ATA' => 'I', # IsoleTcine
			'ATC' => 'I', # IsoleTcine
			'ATT' => 'I', # IsoleTcine
			'ATG' => 'M', # Methionine
			'ACA' => 'T', # Threonine
			'ACC' => 'T', # Threonine
			'ACG' => 'T', # Threonine
			'ACT' => 'T', # Threonine
			'AAC' => 'N', # Asparagine
			'AAT' => 'N', # Asparagine
			'AAA' => 'K', # Lysine
			'AAG' => 'K', # Lysine
			'AGC' => 'S', # Serine
			'AGT' => 'S', # Serine
			'AGA' => 'R', # Arginine
			'AGG' => 'R', # Arginine
			'GTA' => 'V', # Valine
			'GTC' => 'V', # Valine
			'GTG' => 'V', # Valine
			'GTT' => 'V', # Valine
			'GCA' => 'A', # Alanine
			'GCC' => 'A', # Alanine
			'GCG' => 'A', # Alanine
			'GCT' => 'A', # Alanine
			'GAC' => 'D', # Aspartic Acid
			'GAT' => 'D', # Aspartic Acid
			'GAA' => 'E', # GlTtamic Acid
			'GAG' => 'E', # GlTtamic Acid
			'GGA' => 'G', # Glycine
			'GGC' => 'G', # Glycine
			'GGG' => 'G', # Glycine
			'GGT' => 'G'  # Glycine
			);
	##### WE WILL GET ALL HUMAN SEQUENCES #####
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
	
	#print Dumper(%sequences);

	###### First step --> loop through all transcripts found in that type of cancer ###########
	#open (TRANSCRIPT_LIST, "/Users/opich/CRG/Projects/NegSel_Cancer/TCGA_Split/$cancer_type/transcripts.list") or die print "$cancer_type";
	open (TRANSCRIPT_LIST, $ARGV[1]) or die "caca en transcript";
	while (<TRANSCRIPT_LIST>)
	{
		my $read = $_; 
		chomp $read;
	########################################################################################################
	###################### LETS TAKE A LOOK AT CpGs first!!! ###############################################
	########################################################################################################
		#print Dumper($read,$sequences{$read});
	
		my @sequence_trans = split ('', $sequences{$read});
		my %cpg_positions;
		my @normal_positions;
		my %check_cpg;
	##################### LETS DEFINE TWO LIST, ONE FOR THE CpG POSITIONS, THE OTHER FOR THE OTHER POSITIONS
		for my $i (0 .. $#sequence_trans)
		{
	      	if ($sequence_trans[$i] eq "C" && $sequence_trans[$i+1] eq "G")
	      	{
	      		if (exists $cpg_positions{$i+1})
	      		{}
	      		else
	      		{
	      			#push (@cpg_positions, $i+1);
	      			$cpg_positions{$i+1} = 1;
	      		}
	      		if (exists $cpg_positions{$i+2})
	      		{}
	      		else
	      		{
	      			$cpg_positions{$i+2} = 1;
	      		}
	      	}
	      	if (exists $check_cpg{$i+1})
	      	{

	      	}
	      	else
	      	{
	      		push (@normal_positions, $i+1);
	      	}
		}

	############ NOW WE HAVE ALL THE POSITIONS FOR WHICH A CPG IS FOUND AND THE OTHER POSITIONS #############
		#open (OUTFILE, ">$read.new") or die;
		########## NOW LETS OPEN EACH TRANSCRIPT FILE AND GET ALL THE MUTATIONS ############## 
#		open (TRANSCRIPT_MUTATED, "$read") or die;
		my %aloha;
		my %final_check;
#		my %syn_mut; my %syn_mut_only; my %nonsyn_mut; my %nonsyn_mut_only; my %syn_mut_cpg; my %nonsyn_mut_cpg; my %syn_mut_cpg; my %syn_mut_only_cpg; my %nonsyn_mut_only_cpg;
=cut
		while (<TRANSCRIPT_MUTATED>)
		{
			my $llegeix = $_; 
			my @cancer = split ('\s+', $llegeix);
			my $ID = @cancer[0];
			my ($sample) = ($llegeix =~ m/TCGA\-(\S+)/);
			my @sample_id = split ('-', $sample);
			my $tumour_id = join ('_', @sample_id[0], @sample_id[1],@sample_id[2]); #### LETS MAKE SURE NO REPEATED PATIENTS ARE FOUND
			my $uniq_id = join ('_', $ID,$tumour_id);
			my $position_cancer = @cancer[8];
			my ($analysis)=($position_cancer =~ m/(\S+)\/\S+/); ########## LETS GET THE POSITION WHERE THE CHANGE OCCURS!!!!
			my ($ratio) = ($llegeix =~ m/\S+_\S+_(\S+)/); #### TYPE OF MUTATION!!!!!!
			if (exists ($final_check{$uniq_id})) #### MIRAR QUE NOMES HI HAGI UN PER PACIENT!!!!
			{}
			else
			{
				#LETS CHECK IF IT IS FOUND IN A POSITION LABELED AS A CpG#
				if (exists $cpg_positions{$analysis})
				{
					if ($llegeix =~ m/synonymous_variant/ )
					{
						$syn_mut{cpg} = $syn_mut{cpg} +1; ## cada tipus de mutation va al hash i es suma un 
						if (exists $aloha{$ID})
						{
						}
						else
						{
							$syn_mut_only{cpg} = $syn_mut_only{cpg} +1; ## aquest es l'unic!
							$aloha{$ID} = 1;
						}
					}
					elsif ($llegeix=~ m/missense_variant/)
					{
						$nonsyn_mut{cpg} = $nonsyn_mut{cpg} +1; ## cada tipus de mutation va al hash i es suma un 
						if (exists $aloha{$ID})
						{
						}
						else
						{
							$nonsyn_mut_only{cpg} = $nonsyn_mut_only{cpg} +1; ## aquest es l'unic!
							$aloha{$ID} = 1;
						}
					}
					$final_check{$uniq_id} = 1;
				}
				else
				{

					if ($llegeix =~ m/synonymous_variant/ )
					{
						$syn_mut{$ratio} = $syn_mut{$ratio} +1; ## cada tipus de mutation va al hash i es suma un 
						if (exists $aloha{$ID})
						{
						}
						else
						{
							$syn_mut_only{$ratio} = $syn_mut_only{$ratio} +1; ## aquest es l'unic!
							$aloha{$ID} = 1;
						}
					}
					elsif ($llegeix=~ m/missense_variant/)
					{
						$nonsyn_mut{$ratio} = $nonsyn_mut{$ratio} +1; ## cada tipus de mutation va al hash i es suma un 
						if (exists $aloha{$ID})
						{
						}
						else
						{
							$nonsyn_mut_only{$ratio} = $nonsyn_mut_only{$ratio} +1; ## aquest es l'unic!
							$aloha{$ID} = 1;
						}
					}
					$final_check{$uniq_id} = 1;
				}
			}
		}
=cut
		
		my @aa = ($sequences{$read}=~ m/.../g );
		#print Dumper(@aa);
		###### LETS GET ALL THE SITES FOR EACH SPECIFIC VARIATION ####
		my %sites_ratio_syn;
		my %sites_ratio_nonsyn;
		my $syn;
		my $nonsyn;
		my $nucleotide_count = 0;
		foreach my $triplet (@aa)
		{
			
			my @list = split ('', $triplet);
			my $nucl1 = $list[0]; #first element of the triplet
			my $nucl2 = $list[1]; #second element of the triplet
			my $nucl3 = $list[2]; #third element of the triplet
			$nucleotide_count++;
			if (exists $cpg_positions{$nucleotide_count}) ######### THAT MEANS THAT THAT SITE IS A CpG, SO WE SEPARATE IT!!!!!
			{
				my $possible = $variations{$nucl1}; #all possible mutations for the first nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl1, $thing);
					my $new_codon = join ('', $thing, $nucl2, $nucl3);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{cpg} = $sites_ratio_nonsyn{cpg} +1;
					}
					else
					{
						$sites_ratio_syn{cpg} = $sites_ratio_syn{cpg} +1;
					}
				}
			}
			else
			{
				my $possible = $variations{$nucl1}; #all possible mutations for the first nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl1, $thing);
					my $new_codon = join ('', $thing, $nucl2, $nucl3);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{$actual_change} = $sites_ratio_nonsyn{$actual_change} +1;
					}
					else
					{
						$sites_ratio_syn{$actual_change} = $sites_ratio_syn{$actual_change} +1;
					}
				}
			}
			$nucleotide_count++;
			if (exists $cpg_positions{$nucleotide_count})
			{
				my $possible = $variations{$nucl2}; #all possible mutations for the second nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl2, $thing);
					my $new_codon = join ('', $nucl1, $thing, $nucl3);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{cpg} = $sites_ratio_nonsyn{cpg} +1;
					}
					else
					{
						$sites_ratio_syn{cpg} = $sites_ratio_syn{cpg} +1;
					}
				}
			}
			else
			{
				my $possible = $variations{$nucl2}; #all possible mutations for the second nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl2, $thing);
					my $new_codon = join ('', $nucl1, $thing, $nucl3);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{$actual_change} = $sites_ratio_nonsyn{$actual_change} +1;
					}
					else
					{
						$sites_ratio_syn{$actual_change} = $sites_ratio_syn{$actual_change} +1;
					}
				}
			}
			$nucleotide_count++;
			if (exists $cpg_positions{$nucleotide_count})
			{
				my $possible = $variations{$nucl3};  #all possible mutations for the third nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl3, $thing);
					my $new_codon = join ('',  $nucl1, $nucl2, $thing);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{cpg} = $sites_ratio_nonsyn{cpg} +1;
					}
					else
					{
						$sites_ratio_syn{cpg} = $sites_ratio_syn{cpg} +1;
					}
				}	
			}
			else
			{
				my $possible = $variations{$nucl3};  #all possible mutations for the third nucleotide(so A to C,G,T and so on)
				my @all_possible = split ('', $possible);
				foreach my $thing (@all_possible)
				{
					my $actual_change = join ('/', $nucl3, $thing);
					my $new_codon = join ('',  $nucl1, $nucl2, $thing);
					if ($genetic_code{$triplet} ne $genetic_code{$new_codon})
					{
						$sites_ratio_nonsyn{$actual_change} = $sites_ratio_nonsyn{$actual_change} +1;
					}
					else
					{
						$sites_ratio_syn{$actual_change} = $sites_ratio_syn{$actual_change} +1;
					}
				}	
			}
		}		
		my %complement_exists;

		#print "\n$read\t";
			foreach my $key(sort keys %sites_ratio_nonsyn)
			{
				my @transverse = split ('/', $key);
				my $complementary_change = $DNA_complement{$transverse[0]}.'/'.$DNA_complement{$transverse[1]};
				if (exists $complement_exists{$key})
				{}
				else
				{
				#	my $nonsyn_uniq = $nonsyn_mut_only{$key} + $nonsyn_mut_only{$complementary_change};
				#	my $syn_uniq = $syn_mut_only{$key} + $syn_mut_only{$complementary_change};
				#	my $nonsyn_total =  $nonsyn_mut{$key} + $nonsyn_mut{$complementary_change};
				#	my $syn_total =  $syn_mut{$key} + $syn_mut{$complementary_change};
					my $nonsyn_sites_total = ($sites_ratio_nonsyn{$key} + $sites_ratio_nonsyn{$complementary_change})/3;
					my $syn_sites_total = ($sites_ratio_syn{$key} + $sites_ratio_syn{$complementary_change})/3;
					#print OUTFILE"$key\t$complementary_change\t$nonsyn_total\t$syn_total\t$nonsyn_sites_total\t$syn_sites_total\t$nonsyn_uniq\t$syn_uniq\n"; 
					#print PROVA "$key\t$complementary_change\t$nonsyn_total\t$syn_total\t$nonsyn_sites_total\t$syn_sites_total\t$nonsyn_uniq\t$syn_uniq\t";
					if ($key =~ m/A\/C|A\/G|A\/T|C\/A|C\/G|C\/T|cpg/){
						print "$read\t$key\t$complementary_change\t$nonsyn_sites_total\t$syn_sites_total\n"; 
					}else{
						print "$read\t$complementary_change\t$key\t$nonsyn_sites_total\t$syn_sites_total\n"; 
					}				

					#print PROVA "$key\t$complementary_change\t$nonsyn_sites_total\t$syn_sites_total\t"; 
					delete $sites_ratio_nonsyn{$complementary_change};
					$complement_exists{$complementary_change} = 1;
				}	
			}
		if (exists $sites_ratio_nonsyn{cpg})
		{
			
		}
		else
		{
			#print OUTFILE "cpg\trev/\t\\0\t\\0\t\\0\t\\0\t\\0\t\\0\n";
			print  $read."\tcpg\trev/\t0\t0\n";
		}

	}

