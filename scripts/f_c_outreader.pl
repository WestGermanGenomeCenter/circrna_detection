#!/usr/bin/perl -w
# parsing the find_circ output : add coordinates in chr:start-end format, put out only important .bed columns

use strict;
use Getopt::Long qw(GetOptions);

my$infile="in.processed";
my$outfile="out.tsv";
my$alternate_chroms=0;# ignore alternate chroms by default
my$strict=0; # set to 1 to only accept 40x40 alignment quality

GetOptions('infile|i=s' => \$infile,'outfile|o=s'=> \$outfile,'ignore_alt|alt=i' => \$alternate_chroms,'strict|s=i'=>\$strict) or warn "Using default parameters: $0 --from NAME\n";

chomp($infile,$outfile,$alternate_chroms,$strict);

# in
open(IN,$infile)||die "$!";
my@newin = <IN>;

# out
open(ND,">",$outfile)|| die "$!";
print ND "coordinates\tstrand\tsampleid\tunique_counts\tscore\tscore\tRefseqID\n";


foreach my $line (@newin){
      chomp $line;
      my@parts=split(/\t+/,$line);	# split line to find coordinates of gene
      my$chr=$parts[0];
      my$beg=$parts[1];
      my$end=$parts[2];
      #print"$chr:$beg-$end\n";
      my$un="$chr:$beg-$end";		# all coordinades together are the unique id here
      chomp $un;
      my$newline="$un\t$line";
      my@all_things=split(/\t+/,$newline);
      my$ccord=$all_things[0]; # should be chr10:101654702-101656154
      my$long_id=$all_things[4]; # should be auto_circ_004447
      $long_id=~s/run_//;
      ## this string needs some work : remove circ_..
      $long_id =~s/circ\_*.[0-9]{1,20}//ig ;
      # removed circ_8945 for each line
      my$strand=$all_things[6];
      my$uniques=$all_things[7];
      my$bestqa=$all_things[8];
      my$bestqb=$all_things[9];
      my$refseqid=$all_things[20];
      # check tRefseqID
      if(!($refseqid=~/N/)){
        $refseqid=$all_things[21];
        if(!($refseqid=~/N/)){
          $refseqid=$all_things[22]; # gene name as last option
        }
        # then one of the othe two columns can be the refseqid- 23 or 24
      }
      if($strict){
        if($bestqa < 40 || $bestqb < 40 ){ # here we only want 40x40 alignments, else skip
          next;
        }
      }
      if($alternate_chroms =~/1/gx){
            print ND "$ccord\t$strand\t$long_id\t$uniques\t$bestqa\t$bestqb\t$refseqid\n";
      }
      else{
            if(!($ccord=~/chrUn_gl/)){ # check for unsure coordinates
                  print ND "$ccord\t$strand\t$long_id\t$uniques\t$bestqa\t$bestqb\t$refseqid\n";
            }
      }
}
1;
