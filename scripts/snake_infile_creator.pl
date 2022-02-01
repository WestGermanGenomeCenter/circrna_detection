#!/usr/bin/perl -w
use strict;
use Getopt::Long qw(GetOptions);

my$inputfile="fastq_list.txt";
my$group_string="default_group";
my$to_replace="001";
my$replace_into="002";
my$dir="";


GetOptions('to_replace|l1=s' => \$to_replace,'infile|i=s' => \$inputfile,'replace_into|l2=s' => \$replace_into,'dir|d=s'=>\$dir) or warn "Using default parameters: $0 --from NAME\n";

chomp $to_replace;
chomp $replace_into;
chomp $inputfile;
chomp $dir;

open(IN,$inputfile)|| die "$!";	# infile is a .csv file steptwo output.csv
my@lines=<IN>;
my$i=0;
print "samples\n";
foreach my $singleline (@lines){
      my$i++;
      # check for first file
      # if there use parameters to replace into second file
      # make nice samplename,
      # print in stdout: first file, second file, sample name, group
      chomp $singleline;
      my$file_1_check=`ls -1 $dir/$singleline`;
            if(($file_1_check=~/\.fastq/)&&($singleline=~/$to_replace/)){ # if file is there and the current line has what every second file should have
                  # first file is there
                  chomp $singleline;
                  my$to_change=$singleline;
                  my$to_change_2=$singleline;
                  $to_change =~s /$to_replace/ooo/o;# so we cant do it all at once, need some very unique thing as surrogate
                  $to_change=~s/ooo/$replace_into/o;
                  #print "done changing file $singleline into $to_change\n";
                  my$file_two_check=`ls -1 $dir/$to_change`;
                  if($file_two_check=~/\.fastq/){
                        # make samplename nice
                        $to_change_2=~s/\.fastq//;
                        $to_change_2=~s/\.fq//;
                        $to_change_2=~s/\.fastq.gz//;
                        $to_change_2=~s/\.fq.gz//;
                        $to_change_2=~s /$to_replace//;

                        # now we have one thing per line - the sample anem without file endings and without lane identifier
                        # we still only print if both files are there in the first place


                                              print "$to_change_2\n";
                  }

      }
}
1;
