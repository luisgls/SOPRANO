#!/bin/bash
##Luis Zapata 2018. SOPRANO: Calculate selection using dN/dS in target regions
###Important to edit before running
####Hardcode where genome and fasta file are
####Version 3.1

BASEDIR=/home/lortiz/tools/SOPRANO
SUPA=$BASEDIR/data 
TRANS=$BASEDIR/data/ensemble_transcriptID.fasta
TMP=$BASEDIR/tmp/
FASTA=/location/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
GENOME=/location/to/Homo_sapiens.GRCh37.75.dna.primary_assembly.genome


###Check arguments before running
if (($# < 8));  then
	echo "Please provide the four mandatory arguments and optionals (default in parenthesis)
	 -i Annotation - VEP annotated file
         -b Bed file - Provide bed file with Protein coordinates named by Transcript (ENSTXXXXXX 123 135)
         -o Results - Directory
         -n Name - Give name to results
         Optional arguments:
	 -m (ssb7)/ssb192
         -r (target)/random - Calculate a dNdS value for a random region similar to the target
         -e (true)/false - Exclude driver genes from calculation
         -t (bed file with regions to randomize)
"
	exit 1
fi	


##Read all arguments
echo "STEP 1: Reading all arguments and input files"
while getopts "i:b:o:n:m:r:e:t:" opt; do
    case $opt in
	i)
	echo "#-i was triggered, Input File: $OPTARG" >&2 
	FILE=$OPTARG
	;;
        o)
        echo "#-o was triggered, Output folder: $OPTARG" >&2
        OUT=$OPTARG
        ;;
        b)
        echo "#-b was triggered, Target Protein Regions: $OPTARG" >&2
        BED=$OPTARG
        ;;
        n)
        echo "#-n was triggered, Name: $OPTARG" >&2
        NAME=$OPTARG
        ;;
	m)
	echo "#-m was triggered, Substitution model: $OPTARG" >&2 
	MUTRATE=$OPTARG
	;;
        r)
        echo "#-r was triggered, Parameter: $OPTARG" >&2
        MODEL=$OPTARG
        ;;
        e)
        echo "#-e was triggered, Parameter: $OPTARG" >&2
        DRIVER=$OPTARG
        ;;
        t)
        echo "#-t was triggered, Parameter: $OPTARG" >&2
        TARGET=$OPTARG
        ;;
	\?)
	echo "#Invalid option: -$OPTARG" >&2
	exit 1
	;;
	:)
	echo "#Option -$OPTARG requires an argumrent." >&2
	exit 1
	;;	
    esac
done

###Check if file exists and is not empty
if [[ ! -s $FILE ]] ; then
    echo "File $FILE is not there or is empty, aborting."
    exit 1
fi
if [[ ! -s $BED ]] ; then
    echo "File $BED is not there or is empty, aborting."
    exit 1
fi

echo "STEP 2: Preparing coordinate files to intersect data"

###Get list of transcript from annotated Bed file, filter out those transcripts not present in the database.
cut -f1 $BED | sort -u | fgrep -w -f - $SUPA/ensemble_transcript_protein.length > $TMP/$NAME.protein_length_filt.txt
cut -f1 $BED | sort -u | fgrep -w -f - $SUPA/ensemble_transcript.length > $TMP/$NAME.transcript_length_filt.txt

###Randomize protein positions of target region
if [[ $MODEL = "random" ]];
then
    echo "#Option random enabled, running dN/dS for matching target region length"
    ## Define excluded regions for the randomization
    ## Excluding region to be tested
    rm $TMP/$NAME.exclusion.ori
    cut -f1,2,3 $BED > $TMP/$NAME.exclusion.ori
    ## Excluding the first two aminoacids of the transcript to be tested
    cut -f1 $BED | awk '{OFS="\t"}{print $1,0,2}' | sortBed -i stdin  >> $TMP/$NAME.exclusion.ori

    ##Sort excluded regions file
    sortBed -i $TMP/$NAME.exclusion.ori > $TMP/$NAME.exclusion.bed
    bedtools shuffle -i $BED -g $TMP/$NAME.protein_length_filt.txt -excl $TMP/$NAME.exclusion.bed -chrom > $TMP/$NAME.epitopes.ori2

    ###Option to randomize only on a set of target protein regions
    if [ -s "$TARGET" ]
    then
        echo "Target file to randomize regions provided"
        ###Subset transcript protein length and transcript length
        cut -f1 $TARGET | sort -u | fgrep -w -f - $SUPA/ensemble_transcript_protein.length >> $TMP/$NAME.protein_length_filt.txt
        cut -f1 $TARGET | sort -u | fgrep -w -f - $SUPA/ensemble_transcript.length >> $TMP/$NAME.transcript_length_filt.txt
        bedtools shuffle -i $BED -g $TMP/$NAME.protein_length_filt.txt -incl $TARGET -noOverlapping > $TMP/$NAME.epitopes.ori2
    else
        echo "Target file to randomize regions not provided, using default (all)"
    fi
else
    ## If non randomized
    sort -u $BED > $TMP/$NAME.epitopes.ori2
    echo "#Calculating dN/dS for target region"
fi

##Exclude positively selected genes
if [[ $DRIVER = "false" ]];
    then
        echo "#Option excluding driver disabled"
        cp $TMP/$NAME.epitopes.ori2 $TMP/$NAME.epitopes.bed
else
    ##Exclude list of driver transcript genes provided in supplementary data
    fgrep -w -v -f $SUPA/genes2exclude.txt $TMP/$NAME.epitopes.ori2 > $TMP/$NAME.epitopes.bed
    echo "#Option excluding driver enabled (default)"
fi

##Get the complement for the rest of the protein (nonepitope) and make sure we dont look into other proteins not present in the epitope file
sortBed -i $TMP/$NAME.epitopes.bed | complementBed -i stdin -g $TMP/$NAME.protein_length_filt.txt > $TMP/$NAME.intra_epitopes_prot.tmp
cut -f1 $TMP/$NAME.epitopes.bed | sort -u | fgrep -w -f - $TMP/$NAME.intra_epitopes_prot.tmp > $TMP/$NAME.intra_epitopes_prot.bed

## transform protein coordinates to CDS coordinates and get the complement at the transcript level
#Line to get the transcript to be calculated for ssb192
if [[ $MUTRATE = "ssb192" ]];
    then
        awk '{OFS="\t"}{if( (($2*3)-6) >= 0 ) {print $1,($2*3)-6,$3*3+3,$0} else{print $1,($2*3)-3,$3*3+3,$0}}' $TMP/$NAME.epitopes.bed > $TMP/$NAME.epitopes_cds.bed
else
        awk '{OFS="\t"}{print $1,($2*3)-3,$3*3,$0}' $TMP/$NAME.epitopes.bed > $TMP/$NAME.epitopes_cds.bed
fi
    
sortBed -i $TMP/$NAME.epitopes_cds.bed | complementBed -i stdin -g $TMP/$NAME.transcript_length_filt.txt > $TMP/$NAME.intra_epitopes.tmp
cut -f1 $TMP/$NAME.epitopes.bed | sort -u | fgrep -w -f - $TMP/$NAME.intra_epitopes.tmp > $TMP/$NAME.intra_epitopes_cds.bed

##get transcript sequence or epitopes and nonepitopes
## get fasta sequence for each epitope peptide

###
echo "STEP 3: Obtaining target and non-target fasta regions"
bedtools getfasta -fi $TRANS -bed $TMP/$NAME.epitopes_cds.bed -fo $TMP/$NAME.epitopes_cds.fasta
bedtools getfasta -fi $TRANS -bed $TMP/$NAME.intra_epitopes_cds.bed -fo $TMP/$NAME.intra_epitopes_cds.fasta

#List of transcript:regions to estimate number of sites
grep ">" $TMP/$NAME.epitopes_cds.fasta | sed 's/>//g' > $TMP/$NAME.listA
grep ">" $TMP/$NAME.intra_epitopes_cds.fasta | sed 's/>//g' > $TMP/$NAME.listB


echo "STEP 4: Running analysis"
####################################### Running with SSB192 #############################################
if [[ $MUTRATE = "ssb192" ]];
    then
        echo "#Estimate all theoretically possible 192 substitutions in target and non-target regions"
        perl $BASEDIR/scripts/calculate_sites_signaturesLZ_192.pl $TMP/$NAME.epitopes_cds.fasta $TMP/$NAME.listA > $TMP/$NAME.listA.sites
        perl $BASEDIR/scripts/calculate_sites_signaturesLZ_192.pl $TMP/$NAME.intra_epitopes_cds.fasta $TMP/$NAME.listB > $TMP/$NAME.listB.sites
        
        #Sum all possible accross the target region and the non-target region
        awk '{print "test_"$2"_"$3"\t0\t1\t"$0}' $TMP/$NAME.listA.sites | sortBed -i stdin |
        mergeBed -i stdin -c 7,8 -o sum,sum | cut -f1,4,5 > $TMP/$NAME.listA.totalsites
        awk '{print "test_"$2"_"$3"\t0\t1\t"$0}' $TMP/$NAME.listB.sites | sortBed -i stdin |
        mergeBed -i stdin -c 7,8 -o sum,sum | cut -f1,4,5 > $TMP/$NAME.listB.totalsites
        
        ##Get columns to process and obtain trinucleotide context from file
        echo "Proccessing $FILE VEP annotated file to estimate 192 rate parameters..."
        
        
        #1
        echo "In case file starts with column 1 in alternative format"
        perl $BASEDIR/scripts/fixsimulated.pl $FILE > $NAME
        cut -f1,5,7,11,12 $NAME | grep -v "#" | sed 's/_/\t/1' | sed 's/_/\t/1' | awk -F"\t" '{OFS="\t"}{if($6!="-"&&length($3)==3){print $1,$2-1,$2,$0}}' > $NAME.tmp
        #2 For context correction
        bedtools slop -i $NAME.tmp -b 1 -g $GENOME | bedtools getfasta -fi $FASTA -bed stdin -tab | sed 's/:/\t/g' | sed 's/-/\t/g' > $NAME.tmp.bed
        #3
        paste $NAME.tmp $NAME.tmp.bed | cut -f6,14 - |  awk -F "/" '{FS="/"}{OFS="\t"}{print $1,$2}' |  awk -F "" '{FS=""}{OFS="\t"}{if( ($1==$6) && ($3!="-") ){print "GOOD"}else{print "FAIL"}}' > $NAME.flag
        #4
        paste $NAME.tmp $NAME.tmp.bed $NAME.flag | awk '{if($15=="GOOD"){print $0}}' - | cut -f6,14 - |  sort -k2,2 |uniq -c | sed 's/^ \+//g' | sort -k1,1 -n | sed 's/ /\t/g' | awk '{OFS="\t"}{print $3,$2,$1}' | sed -e 's/\t[A-Z]\//_/g' > $NAME.finalVEP.triplets.counts
        
        
        if [ -s "$COUNTS" ] # NOTE: Is this implemented??
        then 
                cat $COUNTS > $NAME.finalVEP.triplets.counts
                echo "Using provided file $COUNTS for correction"
        else
                echo "Using estimated rate parameters for correction"
        fi
        
                echo "Estimating 192 rate parameters"
        #        perl $BASEDIR/scripts/transform192to7.pl $NAME.finalVEP.triplets.counts data/final_translate_SSB192toSSB7 | awk -F "\t" '{OFS="\t"}{print $3,1,2,$2}' | sortBed -i stdin | mergeBed -i stdin -c 4 -o sum | awk '{OFS="\t"}{print "Estimated",$1,$4}' | sed 's/_/\//g' > tmp_to_7
                cp $NAME.finalVEP.triplets.counts $NAME.finalVEP.triplets192.counts
        #	mv tmp_to_7 $NAME.finalVEP.triplets.counts
        
        
        ##Check if triplet counts exist and how many possibilities there were (expected 192 or 7)
        if [ -s "$NAME.finalVEP.triplets.counts" ] 
        then
                muts=`wc -l $NAME.tmp.bed | awk '{ print $1 }'`
                back=`wc -l $NAME.finalVEP.triplets.counts | awk '{ print $1 }'`
                fails=`grep -c "FAIL" $NAME.flag`
                echo "Rate parameter file $NAME.finalVEP.triplets.counts has $back lines of data."
                rm $NAME.tmp $NAME.tmp.bed $NAME.flag 
                echo "Proccesed $muts mutations from VEP file. $fails mutations were discarded (indels or reference conflicts)"
                
            if [ "$back" -lt 7 ]
                then
                        echo "Very few rate parameters to analyse data, please provide rate parameter file with option -c"
                        exit 1
                    fi
                    
                
        else
                echo "$NAME.finalVEP.triplets.counts is empty."
                exit 1
        fi
        
 
        
        
        echo "#Make correction based on total sites" 
        perl $BASEDIR/scripts/correct_update_epitope_sitesV3.pl $TMP/$NAME.listA.totalsites $NAME.finalVEP.triplets.counts > $TMP/$NAME.final_corrected_matrix_A.txt
        perl $BASEDIR/scripts/correct_update_epitope_sitesV3.pl $TMP/$NAME.listB.totalsites $NAME.finalVEP.triplets.counts  > $TMP/$NAME.final_corrected_matrix_B.txt

####################################END of running with SSB192###############################################
else
####################################### Running with SSB7 #############################################
    echo "#Estimate all theoretically possible 7 substitutions in target and non-target regions"
    perl $BASEDIR/scripts/calculate_sites_signaturesLZ.pl $TMP/$NAME.epitopes_cds.fasta $TMP/$NAME.listA > $TMP/$NAME.listA.sites
    perl $BASEDIR/scripts/calculate_sites_signaturesLZ.pl $TMP/$NAME.intra_epitopes_cds.fasta $TMP/$NAME.listB > $TMP/$NAME.listB.sites
    
    #Sum all possible accross the target region and the non-target region
    awk '{print "test_"$2"_"$3"\t0\t1\t"$0}' $TMP/$NAME.listA.sites | sortBed -i stdin |
    mergeBed -i stdin -c 7,8 -o sum,sum | cut -f1,4,5 > $TMP/$NAME.listA.totalsites
    awk '{print "test_"$2"_"$3"\t0\t1\t"$0}' $TMP/$NAME.listB.sites | sortBed -i stdin |
    mergeBed -i stdin -c 7,8 -o sum,sum | cut -f1,4,5 > $TMP/$NAME.listB.totalsites
    
    ##Get columns to process and obtain trinucleotide context from file
    echo "Proccessing $FILE VEP annotated file to estimate 7 rate parameters..."
    #1 Transform the ENST from simulated mutations coming frmo SISMO to chr pos start using eprl script
    
    echo "In case file starts with column 1 in alternative format"
    perl $BASEDIR/scripts/fixsimulated.pl $FILE > $NAME
    
    cut -f1,5,7,11,12 $NAME | grep -v "#" | sed 's/_/\t/1' | sed 's/_/\t/1' | awk -F"\t" '{OFS="\t"}{if($6!="-"&&length($3)==3){print $1,$2-1,$2,$0}}' > $NAME.tmp
    #2 For context correction
    bedtools slop -i $NAME.tmp -b 1 -g $GENOME | bedtools getfasta -fi $FASTA -bed stdin -tab | sed 's/:/\t/g' | sed 's/-/\t/g' > $NAME.tmp.bed
    #3
    paste $NAME.tmp $NAME.tmp.bed | cut -f6,14 - |  awk -F "/" '{FS="/"}{OFS="\t"}{print $1,$2}' |  awk -F "" '{FS=""}{OFS="\t"}{if( ($1==$6) && ($3!="-") ){print "GOOD"}else{print "FAIL"}}' > $NAME.flag
    #4
    paste $NAME.tmp $NAME.tmp.bed $NAME.flag | awk '{if($15=="GOOD"){print $0}}' - | cut -f6,14 - |  sort -k2,2 |uniq -c | sed 's/^ \+//g' | sort -k1,1 -n | sed 's/ /\t/g' | awk '{OFS="\t"}{print $3,$2,$1}' | sed -e 's/\t[A-Z]\//_/g' > $NAME.finalVEP.triplets.counts
    
    
    if [ -s "$COUNTS" ]
    then 
            cat $COUNTS > $NAME.finalVEP.triplets.counts
            echo "Using provided file $COUNTS for correction"
    else
            echo "Using estimated rate parameters for correction"
    fi
    
            echo "Estimating 7 rate parameters"
            perl $BASEDIR/scripts/transform192to7.pl $NAME.finalVEP.triplets.counts $BASEDIR/data/final_translate_SSB192toSSB7 | awk -F "\t" '{OFS="\t"}{print $3,1,2,$2}' | sortBed -i stdin | mergeBed -i stdin -c 4 -o sum | awk '{OFS="\t"}{print "Estimated",$1,$4}' | sed 's/_/\//g' > tmp_to_7
            cp $NAME.finalVEP.triplets.counts $NAME.finalVEP.triplets192.counts
            mv tmp_to_7 $NAME.finalVEP.triplets.counts
    
    
    ##Check if triplet counts exist and how many possibilities there were (expected 192 or 7)
    if [ -s "$NAME.finalVEP.triplets.counts" ] 
    then
            muts=`wc -l $NAME.tmp.bed | awk '{ print $1 }'`
            back=`wc -l $NAME.finalVEP.triplets.counts | awk '{ print $1 }'`
            fails=`grep -c "FAIL" $NAME.flag`
            echo "Rate parameter file $NAME.finalVEP.triplets.counts has $back lines of data."
            rm $NAME.tmp $NAME.tmp.bed $NAME.flag 
            echo "Proccesed $muts mutations from VEP file. $fails mutations were discarded (indels or reference conflicts)"
            
            if [ "$back" -lt 7 ]
                then
                    echo "Very few rate parameters to analyse data, please provide rate parameter file with option -c"
                    exit 1
            fi
    else
            echo "$NAME.finalVEP.triplets.counts is empty."
            exit 1
    fi
    
    echo "#Make correction based on total sites" 
    perl $BASEDIR/scripts/correct_update_epitope_sitesV2.pl $TMP/$NAME.listA.totalsites $NAME.finalVEP.triplets.counts > $TMP/$NAME.final_corrected_matrix_A.txt
    perl $BASEDIR/scripts/correct_update_epitope_sitesV2.pl $TMP/$NAME.listB.totalsites $NAME.finalVEP.triplets.counts  > $TMP/$NAME.final_corrected_matrix_B.txt

fi
####################################### Finish running with SSB7 #############################################

echo "STEP 5: Intersecting information from variants and target regions"
##By frequency of mutations
awk '{NA+=$4}{NS+=$6}END{print NA"\t"NS}' $TMP/$NAME.final_corrected_matrix_A.txt > $TMP/$NAME.epitope_NaNs.txt
awk '{NA+=$4}{NS+=$6}END{print NA"\t"NS}' $TMP/$NAME.final_corrected_matrix_B.txt > $TMP/$NAME.nonepitope_NaNs.txt


#####Get variant counts from vep annotated file, split into silent and nonsilent
##silent
egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|non_coding|splice_acceptor_variant|splice_donor_variant|TF_binding_site_variant|upstream|incomplete|regulatory_region_variant|retained|\?' $FILE | grep -w synonymous_variant |
awk '{if(length($3)>1||$10=="-"){}else{print}}' | cut -f4,5,7,10,89 -  | sed 's/\//\t/g' | awk '{print $2"\t"$4"\t"$4"\t"$3}' |  egrep -v -w -e "coding_sequence_variant" | grep -v "ENSEMBLTRANSCRIPT" > $TMP/$NAME.silent.bed

##nonsilent
egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|non_coding|splice_acceptor_variant|splice_donor_variant|TF_binding_site_variant|upstream|incomplete|regulatory_region_variant|retained|\?' $FILE | grep -w -v synonymous_variant |
awk '{if(length($3)>1||$10=="-"){}else{print}}' | cut -f4,5,7,10,89 -  | sed 's/\//\t/g' | awk '{print $2"\t"$4"\t"$4"\t"$3}' |  egrep -v -w -e "coding_sequence_variant" | grep -v "ENSEMBLTRANSCRIPT" > $TMP/$NAME.nonsilent.bed

##missense only
egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|non_coding|splice_acceptor_variant|splice_donor_variant|TF_binding_site_variant|upstream|incomplete|regulatory_region_variant|retained|\?' $FILE | grep -w -v synonymous_variant | grep -w missense_variant |
awk '{if(length($3)>1||$10=="-"){}else{print}}' | cut -f4,5,7,10,89 -  | sed 's/\//\t/g' | awk '{print $2"\t"$4"\t"$4"\t"$3}' |  egrep -v -w -e "coding_sequence_variant" | grep -v "ENSEMBLTRANSCRIPT" > $TMP/$NAME.missense.bed

##intronic
grep -v "^#" $NAME | grep -w "intron_variant" | grep -v "splice" | awk -F"\t|_" '{FS="\t|_"}{print $1"_"$7"\t"$2"\t"$2"\t"$3}' > $TMP/$NAME.intronic.bed

##Silent
if [ -s "$TMP/$NAME.silent.bed" ]
then
        sils=`wc -l $TMP/$NAME.silent.bed | awk '{ print $1 }'`
        #echo "$sils number of silent mutations in file"
        
else
        echo "file silent mutations empty"
fi

##Nonsilent (nonsense + missense)
if [ -s "$TMP/$NAME.nonsilent.bed" ]
then
        nonsils=`wc -l $TMP/$NAME.nonsilent.bed | awk '{ print $1 }'`
        #echo "$nonsils number of nonsilent mutations in file"
        
else
        echo "file nonsilent mutations empty"
fi

##Missense only
if [ -s "$TMP/$NAME.missense.bed" ]
then
        missense=`wc -l $TMP/$NAME.missense.bed | awk '{ print $1 }'`
        #echo "$missense number of missense mutations in file"
        
else
        echo "file missense mutations empty"
fi

###ntersect different regions from the protein to calculate dNdS
##Check if counts has value larger than 0

if [ -s "$NAME.finalVEP.triplets.counts" ]
then
	echo "Removing rate parameter files"
        rm $NAME.finalVEP.triplets.counts $NAME.finalVEP.triplets192.counts $NAME
else
        echo "WARNING:triplets counts not present"
fi

###Count number of mutations intersecting ON and OFF regions for selected transcripts.

    innonsil=`intersectBed -b $TMP/$NAME.nonsilent.bed -a $TMP/$NAME.epitopes.bed -wo | wc -l | awk '{ print $1 }'`
    inmissen=`intersectBed -b $TMP/$NAME.missense.bed -a $TMP/$NAME.epitopes.bed -wo | wc -l | awk '{ print $1 }'`
    insil=`intersectBed -b $TMP/$NAME.silent.bed -a $TMP/$NAME.epitopes.bed -wo | wc -l | awk '{ print $1 }'`
    outnonsil=`intersectBed -b $TMP/$NAME.nonsilent.bed -a $TMP/$NAME.intra_epitopes_prot.bed -wo | wc -l | awk '{ print $1 }'`
    outmissen=`intersectBed -b $TMP/$NAME.missense.bed -a $TMP/$NAME.intra_epitopes_prot.bed -wo | wc -l | awk '{ print $1 }'`
    outsil=`intersectBed -b $TMP/$NAME.silent.bed -a $TMP/$NAME.intra_epitopes_prot.bed -wo | wc -l | awk '{ print $1 }'`

echo "From a total of $nonsils non-silent, $missense missense-only, and $sils silent mutations in input file:"    
echo "There are $innonsil non-silent, $inmissen missense-only, and $insil silent in the (ON) target region"
echo "There are $outnonsil non-silent, $outmissen missense-only, and $outsil silent in the (OFF) non-target region"

###Modified from nonsilent to missense to calculate for missense only

if [ -s "$TMP/$NAME.data_epitopes" ]
then
        rm $TMP/$NAME.data_epitopes
        echo "INFO: past epitope-associated file, removing"
else
        echo "INFO: intersected file (data_epitopes) not present, creating"
fi

if [ "$insil" -gt 0 ];
then
intersectBed -b $TMP/$NAME.silent.bed -a $TMP/$NAME.epitopes.bed -wo | awk '{OFS="\t"}{print $1,"1","2",$0}' | cut -f1-11 -| sortBed -i stdin | mergeBed -i stdin -c 11 -o count | cut -f1,4 | awk '{print $0"\textra_synonymous_variant"}' >> $TMP/$NAME.data_epitopes
fi

if [ "$innonsil" -gt 0 ];
then
intersectBed -b $TMP/$NAME.nonsilent.bed -a $TMP/$NAME.epitopes.bed -wo | awk '{OFS="\t"}{print $1,"1","2",$0}' | cut -f1-11 -| sortBed -i stdin | mergeBed -i stdin -c 11 -o count | cut -f1,4 | awk '{print $0"\textra_missense_variant"}' >> $TMP/$NAME.data_epitopes
fi

if [ "$outsil" -gt 0 ];
then
intersectBed -b $TMP/$NAME.silent.bed -a $TMP/$NAME.intra_epitopes_prot.bed -wo | awk '{OFS="\t"}{print $1,"1","2",$0}' | cut -f1-11 -| sortBed -i stdin | mergeBed -i stdin -c 10 -o count | cut -f1,4 | awk '{print $0"\tintra_synonymous_variant"}' >> $TMP/$NAME.data_epitopes
fi

if [ "$outnonsil" -gt 0 ];
then
intersectBed -b $TMP/$NAME.nonsilent.bed -a $TMP/$NAME.intra_epitopes_prot.bed -wo | awk '{OFS="\t"}{print $1,"1","2",$0}' | cut -f1-11 - | sortBed -i stdin | mergeBed -i stdin -c 10 -o count | cut -f1,4 | awk '{print $0"\tintra_missense_variant"}' >> $TMP/$NAME.data_epitopes
fi


if [  "$innonsil" -eq 0 ] && [ "$inmissen" -eq 0 ] && [ "$insil" -eq 0 ];
then
    echo "$FILE FAILED: 0 Mutations in target region"
    exit 1
else   
    ### For intronic
    echo "STEP 6:  Calculating dN/dS on target and off target regions"
    intersectBed -a $BASEDIR/data/transcript_intron_length.bed -b $TMP/$NAME.intronic.bed -wo | mergeBed -i stdin -c 4,5,6,10,11 -o mode,mode,mode,collapse,count | awk '{print $4"\t"$8/($6+1)"\t"$8"\t"$6}' >  $TMP/$NAME.intronic.rate
       
    
     #check fot existence of outdir of not create
if [ ! -d $OUT ]
then
	mkdir $OUT
fi

if [ -s $TMP/$NAME.intronic.rate ]
then
    Rscript $BASEDIR/scripts/calculateKaKsEpiCorrected_CI_intron_V3.R $TMP/$NAME.data_epitopes $TMP/$NAME.epitope_NaNs.txt $TMP/$NAME.nonepitope_NaNs.txt $TMP/$NAME.intronic.rate > $OUT/$NAME.SSB_dNdS.txt
else
    Rscript $BASEDIR/scripts/calculateKaKsEpiCorrected_CI.R $TMP/$NAME.data_epitopes $TMP/$NAME.epitope_NaNs.txt $TMP/$NAME.nonepitope_NaNs.txt > $OUT/$NAME.SSB_dNdS.txt
fi

    if [ -s "$OUT/$NAME.SSB_dNdS.txt" ]
    then
            echo "SOPRANO SUCCESS ... removing tmp files"
            rm "$TMP/$NAME."*
            
            echo
    else
            echo "SOPRANO FAILED "
            echo
    fi
fi  
