#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0 --bam bam_file -Q 10 (DEFAULT 0) -q 20 (DEFAULT 10) ";

# declare the options upfront :
my $BAM;


#get options :
GetOptions ("bam=s" => \$BAM    # the mpileup file from the first in pair read
	        	           
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$file_R1 || !$file_R2 || !$out || !$generic) {die $error_sentence}
#=================================

#do not override a file ==========
if(-e $out) { die "File $out Exists, please remove old file or rename output file (--out)"};
#================================= 




my $command_R1 = "samtools view -b -f 64 -b $file > $bam_R1";    
my $command_R2 = "samtools view -b -f 128 -b $file > $bam_R2";

    #run mpileup on the separated Bam files.
my $mpileupR1 = $generic."_R1.mpileup";
my $command_mplipeup_R1 = "samtools mpileup -O -s -q 10 -Q 0 -f $genome $bam_R1 > $mpileupR1";
    
my $mpileupR2 = $generic."_R2.mpileup";
my $command_mplipeup_R2 = "samtools mpileup -O -s -q 10 -Q 0 -f $genome $bam_R2 > $mpileupR2";
    
