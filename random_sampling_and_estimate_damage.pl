#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);
use File::Basename;
my $dirname = dirname(__FILE__);



my $error_sentence = "USAGE : perl $0 --mpileup1 mileup_file_from_first_in_pair_read --mpileup2 mileup_file_from_second_in_pair_read --out damage_file_for_R --id library_name \nOPTIONAL :\n --qualityscore 35 (DEFAULT 30)\n --min_coverage_limit 10 (DEFAULT 1)\n --max_coverage_limit 500 (DEFAULT 100) --soft_masked 1 (DEFAULT 0) -n 100 (DEFAULT 20)\n";

# declare the options upfront :                                                                                                             
my $file_R1;
my $file_R2;
my $out;
my $generic;
my $QUALITY_CUTOFF = 30;
my $COV_MIN = 1;
my $COV_MAX = 100;
my $SOFT_MASKED=0;
my $N=20;
my $SAMPLE = 2000000; #2 Millions positions
#get options :                                                                                                                              
GetOptions ("mpileup1=s" => \$file_R1,    # the mpileup file from the first in pair read                                                  
	    "mpileup2=s" => \$file_R2, #the mpileup file from the second in pair read                                                  
	    "qualityscore=s" => \$QUALITY_CUTOFF,#base quality (either direct from Illumina or recalibrated)                           
	    "out=s" => \$out,#output file name.                                                                                        
	    "id=s" => \$generic,
	    "min_coverage_limit=s" => \$COV_MIN,
	    "max_coverage_limit=s" => \$COV_MAX,
	    "soft_masked=s" => \$SOFT_MASKED,
	    "n=s" => \$N,
	    "sample_size=s" => \$SAMPLE
    ) or die $error_sentence;

#=================================                                                                                                          
#if something went wrong in getting the option, notify :                                                                                    
if (!$file_R1 || !$file_R2 || !$out || !$generic) {die $error_sentence}
#=================================                                                                                                          

#do not override a file ==========                                                                                                          
if(-e $out) { die "File $out Exists, please remove old file or rename output file (--out)"};
#=================================         

for (my $i =0; $i <$N; $i++)
{
    #create tmp files
    my $tmp_R1 = new File::Temp( UNLINK => 1 );
    my $tmp_R2 = new File::Temp( UNLINK => 1 );                                                                                            
    my $c1 = "$dirname/randomized2 $SAMPLE $file_R1 > $tmp_R1";
    my $c2 = "$dirname/randomized2 $SAMPLE $file_R2 > $tmp_R2";
    my $c3 = "perl $dirname/estimate_damage.pl --mpileup1 $tmp_R1 --mpileup2 $tmp_R2 --id $generic --qualityscore $QUALITY_CUTOFF --min_coverage_limit $COV_MIN --max_coverage_limit $COV_MAX --soft_masked $SOFT_MASKED >> $out";
    system($c1);
    system($c2);
    system($c3);
}




