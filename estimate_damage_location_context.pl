#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0 --mpileup1 mileup_file_from_first_in_pair_read --mpileup2 mileup_file_from_second_in_pair_read --out damage_file_for_R --id library_name \nOPTIONAL :\n --qualityscore 35 (DEFAULT 30)\n --min_coverage_limit 10 (DEFAULT 1)\n --max_coverage_limit 500 (DEFAULT 100) --soft_masked 1 (DEFAULT 1) --context 1 (DEFAULT 3)\n";

# declare the options upfront :
my $file_R1;
my $file_R2;
my $out;
my $generic;
my $QUALITY_CUTOFF = 30;
my $COV_MIN = 1;
my $COV_MAX = 100;
my $SOFT_MASKED=1;
my $CONTEXT =3;
#get options :
GetOptions (    "mpileup1=s" => \$file_R1,    # the mpileup file from the first in pair read
	        "mpileup2=s" => \$file_R2, #the mpileup file from the second in pair read
	        "qualityscore=s" => \$QUALITY_CUTOFF,#base quality (either direct from Illumina or recalibrated)
	        "out=s" => \$out,#output file name.
	        "id=s" => \$generic,
	        "min_coverage_limit=s" => \$COV_MIN,
	        "max_coverage_limit=s" => \$COV_MAX,
	        "soft_masked=s" => \$SOFT_MASKED,
		"context=s" => \$CONTEXT
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$file_R1 || !$file_R2 || !$out || !$generic) {die $error_sentence}
#=================================

#need to work on this !!!!!!!!!======================
my $relative_count_R1 = get_relative_count($file_R1);
my $relative_count_R2 = get_relative_count($file_R2);


#open the output file ====================================
open(OUT, ">$out") or die "can't save result into $out\n";
#=========================================================  



foreach my $type (keys %$relative_count_R1)
{
    my $ref = $$relative_count_R1{$type};
    
    foreach my $position (sort {$a<=>$b} keys %$ref)
    {
	my $ref2 = $$ref{$position};
	foreach my $seq (keys %$ref2)
	{
	    
	    my $value1 = $$relative_count_R1{$type}{$position}{$seq}{"relative"};
	    my $value2 = $$relative_count_R2{$type}{$position}{$seq}{"relative"};
	    my $abs1 = $$relative_count_R1{$type}{$position}{$seq}{"absolute"};
	    my $abs2 = $$relative_count_R2{$type}{$position}{$seq}{"absolute"};
	    if ($$relative_count_R1{$type}{$position}{$seq}{"relative"} && $$relative_count_R2{$type}{$position}{$seq}{"relative"})
	    {
		print OUT "$generic\t$type\tR1\t$value1\t$position\t$seq\t$abs1\n";
		print OUT "$generic\t$type\tR2\t$value2\t$position\t$seq\t$abs2\n";
	    }
	}
    }
}

close OUT;

sub clean_bases{
    #remove unwanted sign (insertion / deletion (but retain the + and - and remove the first and last base information)
    my ($pattern, $n)=@_; #for example : ,$,,
    #remove the ^ and the next character that correspond to the first base of the read. 
    $pattern =~ s/\^.//g;
    #remove the dollars sign from the last base of the read. 
    $pattern =~ s/\$//g;
    
    if ($pattern =~/\+(\d+)/ || $pattern=~/\-(\d+)/)
    {
	
	while ($pattern=~/\+(\d+)/)
	{
	    my $size = $1;
	    my $string_to_remove;
	    for(my $i=0; $i<$size; $i++)
	    {
		$string_to_remove = $string_to_remove.".";
	    }
	    $string_to_remove = '.\+\d+'.$string_to_remove;
	    
	    
	    $pattern =~ s/$string_to_remove/\+/;
	    
	}
	while ($pattern=~/\-(\d+)/)
	{
	    my $size = $1;
	    my $string_to_remove;
	    for(my $i=0; $i<$size; $i++)
	    {
		$string_to_remove = $string_to_remove.".";
	    }
	    $string_to_remove = '.\-\d+'.$string_to_remove;
	    
	    
	    $pattern =~ s/$string_to_remove/\-/;
	}	
    }
    return $pattern;
}


sub get_nt{
    my ($line) =@_;
    my ($chr,$loc,$ref, $number, $bases,$q1, $q2, $pos) = split /\t/, $line;
    return ($ref);
}

sub get_relative_count {
    my ($file)=@_; 
    my %result;
    my %final;
    my $l1;
    my $l2;
    my $l3;

    open (FILE, $file) or die;
    while(my $line = <FILE>)
    {
	$line =~ s/\r|\n//g;

	$l1 = $l2;
	$l2 = $l3;
	$l3 = $line;
	if ($l1 && $l2 && $l3){
	    my $nt1 = get_nt($l1);
	    my $nt2 = get_nt($l2);
	    my $nt3 = get_nt($l3);
	    
	    my $seq_all = "no_context";
	    if ($CONTEXT == 1){$seq_all = $nt1."_base"; }
	    if ($CONTEXT == 2){$seq_all = "base_".$nt3; }
	    if ($CONTEXT == 3){$seq_all = $nt1."_base_".$nt3; }

	    my $seq = uc($seq_all);
	    my ($chr,$loc,$ref, $number, $dbases,$q1, $q2, $pos) = split /\t/, $l2;
	    my $bases = clean_bases($dbases, $number); 
	    
	    my @qualities =split //, $q1;	    
	    my @tmp = split //, $bases;
	    my @poss = split/,/, $pos;
	    my $length_mutation = @tmp;

	    if ($SOFT_MASKED ==0)
	    {
		#if we do not want to solf masked the repeat : all reference becomes UPPER CASE. 
		$ref = uc($ref);
	    }
	    if ($number == $length_mutation && $number >= $COV_MIN && $number <= $COV_MAX && $ref =~/[ACTG]/)
	    {
		
		my @base=split (//,$bases);
		for(my $i=0;$i<@base;$i++){

#given a character $q, the corresponding Phred quality can be calculated with:                                             
		    #$Q = ord($q) - 33;  
		    my $quality_score = ord($qualities[$i]) - 33;
		    
		    if ($quality_score > $QUALITY_CUTOFF)
		    {
			my $ch=$base[$i];
			my $rp = $poss[$i];
			if($ch=~/[ATCG]/){
			    my $type = $ref."_".$ch;
			    $result{"type"}{$rp}{$type}{$seq}++;
			    $result{"nt"}{$rp}{$ref}{$seq}++;
			    
			}
			
			elsif($ch eq "."){
			    my $type = $ref."_".$ref;
			    $result{"type"}{$rp}{$type}{$seq}++;
			    $result{"nt"}{$rp}{$ref}{$seq}++;
			    
			    
			}
		    }
		    
		}#end the loop
	    }
	}
    }
    close FILE;
    my $positions = $result{"type"};
    foreach my $pos (keys %$positions)
    {
	my $mutations = $$positions{$pos};
	foreach my $type (keys %$mutations)
	{
	    my $seqs = $$mutations{$type};
	    foreach my $seq (keys %$seqs)
	    {
		my $total_type = $result{"type"}{$pos}{$type}{$seq};
		$type =~/(.)\_(.)/;
		my $total = $result{"nt"}{$pos}{$1}{$seq};
		
		my $value = $total_type/$total;
		
		$final{$type}{$pos}{$seq}{"relative"}=$value;
		$final{$type}{$pos}{$seq}{"absolute"}=$total_type;
	    }
	}
    }
 
    return(\%final);
}


sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCATGCA/;
        return $revcomp;
}
