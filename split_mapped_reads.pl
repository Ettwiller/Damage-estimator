#!/usr/bin/perl
use strict;
use Getopt::Long qw(GetOptions);
use File::Temp qw(tempfile);

my $error_sentence = "USAGE : perl $0 --bam bam_file -genome genome_file.fasta -mpileup1 output_file1 -mpileup2 output_file2 -Q 10 (DEFAULT 0) -q 20 (DEFAULT 10)";

# declare the options upfront :
my $BAM;
my $Q = 0;
my $q = 10;
my $mpileupR1;
my $mpileupR2;
my $genome;
#get options :
GetOptions ("bam=s" => \$BAM,    # the mpileup file from the first in pair read
	    "Q=s" => \$Q,
	    "q=s" => \$q,
	    "mpileup1=s" => \$mpileupR1,
	    "mpileup2=s" => \$mpileupR2,
	    "genome=s" => \$genome
    ) or die $error_sentence;

#=================================
#if something went wrong in getting the option, notify :
if (!$BAM || !$mpileupR1 || !$mpileupR2 || !$genome) {die $error_sentence}
#=================================

#create tmp files for bam1 and bam2 ===================
my $bam_R1 = new File::Temp( UNLINK => 1 );
my $bam_R2 = new File::Temp( UNLINK => 1 );

#do not override a file ==========
#if(-e $out) { die "File $out Exists, please remove old file or rename output file (--out)"};
#================================= 

my $command_R1 = "samtools view -b -f 64 -b $BAM > $bam_R1";    
my $command_R2 = "samtools view -b -f 128 -b $BAM > $bam_R2";

    #run mpileup on the separated Bam files.

my $command_mplipeup_R1 = "samtools mpileup -O -s -q $q -Q $Q -f $genome $bam_R1 > $mpileupR1";
my $command_mplipeup_R2 = "samtools mpileup -O -s -q $q -Q $Q -f $genome $bam_R2 > $mpileupR2";

system($command_R1);
system($command_R2);
system($command_mplipeup_R1);
system($command_mplipeup_R2);

