#!/usr/bin/perl
use strict;






#CONSTANCES
my $QUALITY_CUTOFF = 30;


#THIS PROGRAM IS USED BY THE MAIN PROGRAM DAMAGE-ESTIMATOR
my $file_R1 = $ARGV[0]; #read 1 mpileup file. 
my $file_R2 = $ARGV[1]; #read 2 mpileup file. 
my $generic = $ARGV[2]; #name of the result files containing the damages

if (!$file_R1 || !$file_R2 || !$file_R2) {die "you need to provide read1 mpileup file read2 mpileup file and name of the result files containing the damages\n"};

my %RC;
$RC{'A_C'} = 'T_G';
$RC{'A_A'} = 'T_T';
$RC{'A_T'} = 'T_A';
$RC{'A_G'} = 'T_C';

$RC{'T_C'} = 'A_G';
$RC{'T_A'} = 'A_T';
$RC{'T_T'} = 'A_A';
$RC{'T_G'} = 'A_C';

$RC{'C_C'} = 'G_G';
$RC{'C_A'} = 'G_T';
$RC{'C_T'} = 'G_A';
$RC{'C_G'} = 'G_C';

$RC{'G_C'} = 'C_G';
$RC{'G_A'} = 'C_T';
$RC{'G_T'} = 'C_A';
$RC{'G_G'} = 'C_C';

$RC{'A_-'} = 'T_-';
$RC{'T_-'} = 'A_-';
$RC{'C_-'} = 'G_-';
$RC{'G_-'} = 'C_-';

$RC{'A_+'} = 'T_+';
$RC{'T_+'} = 'A_+';
$RC{'C_+'} = 'G_+';
$RC{'G_+'} = 'C_+';


my $relative_count_R1 = get_relative_count($file_R1);
my $relative_count_R2 = get_relative_count($file_R2);

foreach my $type (keys %RC)
#foreach my $type (keys %$relative_count_R1)
{
    my $ref = $$relative_count_R1{$type};
    foreach my $position (sort {$a<=>$b} keys %$ref)
    {
	my $value1 = $$relative_count_R1{$type}{$position}{"relative"};
	my $value2 = $$relative_count_R2{$type}{$position}{"relative"};
	my $abs1 = $$relative_count_R1{$type}{$position}{"absolute"};
	my $abs2 = $$relative_count_R2{$type}{$position}{"absolute"};
	if ($$relative_count_R1{$type}{$position}{"relative"} && $$relative_count_R2{$type}{$position}{"relative"})
	{
	   

	    print "$generic\t$type\tR1\t$value1\t$position\n";
	    print "$generic\t$type\tR2\t$value2\t$position\n";
	}
    }
}

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





sub get_relative_count {
    my ($file)=@_; 
    my %result;
    my %final;

   
    open (FILE, $file) or die;
    
    while (<FILE>){
	s/\r|\n//g;
	my ($chr,$loc,$ref, $number, $dbases,$q1, $q2, $pos) = split /\t/;
	my $bases = clean_bases($dbases, $number);

	my @tmp = split //, $bases;
	my @poss = split/,/, $pos;
	my $length_mutation = @tmp;
	if ($number == $length_mutation && $number > 0 && $number < 10 && $ref =~/[ACTG]/)
	    #$number == $length_mutation : making sure that the position array and base array are in agreement. 
	    # $number > 0  : position with at least one read
	    #$number < 10  : at 2 M reads across the human genome we are not expecting tp have a high coverage. #!!! be careful of other species.
	    #ref : making sure the position is not soft masked

	{
	    
	    my @base=split (//,$bases);
	    my @qualities =split(//, $q1);
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
			$result{"type"}{$rp}{$type}++;
			$result{"nt"}{$rp}{$ref}++;
			
		    }
		    elsif($ch=~/[\+\-]/ && $rp > 5){ #this is to remove the artefact of alignment algorithm that add indels (?) at the begining. Should be delt differently. 
                        my $type = $ref."_".$ch;
                        $result{"type"}{$rp}{$type}++;
                        $result{"nt"}{$rp}{$ref}++;

                    }

		    elsif($ch eq "."){
			my $type = $ref."_".$ref;
			$result{"type"}{$rp}{$type}++;
			$result{"nt"}{$rp}{$ref}++;
			
		    }
		}
		
	    }#end the loop
	}
    }
    close FILE;
    my $positions = $result{"type"};
    foreach my $pos (keys %$positions)
    {
	my $mutations = $$positions{$pos};
	foreach my $type (keys %$mutations)
	{
	    my $total_type = $result{"type"}{$pos}{$type};
	    $type =~/(.)\_(.)/;
	    my $total = $result{"nt"}{$pos}{$1};
	    
	    my $value = $total_type/$total;
	    
	    $final{$type}{$pos}{"relative"}=$value;
	    $final{$type}{$pos}{"absolute"}=$total_type;
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
