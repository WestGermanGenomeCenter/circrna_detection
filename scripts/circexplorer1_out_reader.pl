#/usr/bin/perl -w
use strict;
# after circexplorer1_starter_1.pl , formats the output of circexplorer into for matrixmaker.pl readable format
# remember that the score will be always the same for each line, as circexplorer gives one score per candidate, and find_circ gives two.
# starting vars
my$currdir=`pwd`;
my$starttime= time;

my$infile=$ARGV[0];
chomp $infile;
my$outfile=$ARGV[1];
chomp $outfile;
my$samplename=$ARGV[2];
chomp $samplename;

open(IN,$infile)|| die "$!";
my@alllines=<IN>;



open(OU,">",$outfile)|| die "$!";

print OU "coordinates\tstrand\tsampleid\tunique_counts\tscore\tscore\tRefseqID\n";
foreach my$single_line (@alllines){
  # split into parts..
  my@lineparts=split(/\t/,$single_line);
  my$chrom=$lineparts[0];
  my$start=$lineparts[1];
  my$end=$lineparts[2];
  my$score=$lineparts[4];
  my$fullcoord="$chrom:$start-$end";

  my$strand=$lineparts[5];
  my$unique_reads=$lineparts[12];
  my$gene_name=$lineparts[14];
  my$read_name=$lineparts[15];
  if($unique_reads >= 2){       # filter parameter for find_circ and DCC: only show circRNA candidates with 2 or more reads
      if(!($chrom=~/chrUn_gl/)){ # check for unsure coordinates 
            my$csvfile="$fullcoord\t$strand\t$samplename\t$unique_reads\t$score\t$score\t$read_name";
            $csvfile=~s/\t\t/\t/g;
            print OU "$csvfile\n";
      }
  }
}
