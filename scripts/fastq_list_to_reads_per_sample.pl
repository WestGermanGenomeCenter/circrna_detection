#!/usr/bin/perl -w
use strict;

use Getopt::Long qw(GetOptions);


my$linfile= "all_fastqs.txt";
my$mode="se";
my$to_replace="001";


GetOptions('to_replace|l1=s' => \$to_replace,'infile|i=s' => \$linfile) or warn "Using default parameters: $0 --from NAME\n";

chomp $linfile;
chomp $mode;
chomp $to_replace;

# go thorough fastq list, if file with lane1 ident is found, wc-l and /4, then save into file



#print "reading input file $linfile ...\n";
open(IN,$linfile)|| die "$!";
my@allelines= <IN>;
print "sample_short\treads_total\n";
foreach my $two_file_line (@allelines){

  if($mode eq "se"){
      chomp $two_file_line;
      if($two_file_line=~/$to_replace/){
        my$reads_file1=`cat $two_file_line|wc -l`;
        chomp $reads_file1;
        my$total_reads_sample= $reads_file1/4;
        my$sample_name=$two_file_line;
        $sample_name=~s /$to_replace//g;# remove lane1/2 indicator
        $sample_name=~s/\.fastq//;# remove file format
      # sample name as tidy as possible
        print "$sample_name\t$total_reads_sample\n";
      }
  }
}
1;
