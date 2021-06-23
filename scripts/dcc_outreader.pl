#/usr/bin/perl -w
use strict;
# starting vars
my$currdir=`pwd`;
my$starttime= time;
# startzed like this ; perl $base_dir/automate_DCC/dcc_outreader.pl CircRNACount_annotated.tsv CircCoordinates processed_run_$samplename.tsv $samplename`;
my$infile=$ARGV[0];# CircRNACount CircRNACount_annotated.tsv file
chomp $infile;
my$linfile=$ARGV[1];# CircCoordinates file
chomp $linfile;
# outfile
my$outfile=$ARGV[2];
chomp $outfile;
my$samplename=$ARGV[3];
chomp $samplename;

open(IN,$infile)|| die "$!";
my@alllines=<IN>;# CircRNACount CircRNACount_annotated.tsv file
open(IZ,$linfile)|| die "$!";
my@allllines=<IZ>;# CircCoordinates file
my%mapping_hash;
foreach my $line (@allllines){# CircCoordinates file
  # get scores and coordinates, match in later file
  my@loineparts=split(/\t/,$line);
  my$score = $loineparts[4];
  my$chrom=$loineparts[0];
  my$start=$loineparts[1];
  my$end=$loineparts[2];
 # my$strand=$loineparts[5];
  my$fullcoordsi="$chrom:$start-$end";
  $mapping_hash{"$fullcoordsi"}="$score";# coords as key, score as value
}
my@allmapco= keys(%mapping_hash);
open(OU,">",$outfile)|| die "$!";
# order has to be the same in those two files, otherwise this will not work
print OU "coordinates\tstrand\tsampleid\tunique_counts\tscore\tscore\tRefseqID\n";
# foreach line i.e circ candidate from dcc
for (my $var = 0; $var < scalar(@alllines); $var++) { #CircRNACount
  my$countline=$alllines[$var];# CircRNACount CircRNACount_annotated.tsv file
  # in the annotated file, we see
  # chr20	25656368	25656448	+	1	chr20	25654850	25677469	NM_015655	0	-	25655667	25667053	0	5	2823,96,127,76,75,	0,11352,11781,12176,22544,
  #DCCcirc coordinates not yet -1 corrected| quant | bedtools found gene coordinates surrounding circRNA coords | refseqid connected to those | exons

  if(!($countline=~/Start/)){ # header check
      my$cordline=$allllines[$var];
      my@count_line_parts=split(/\t/,$countline);# CircRNACount file
      # example line : chr1	782719	782927	-	1	chr1	762970	794826	NR_015368	0	+	794826	794826	0	8	185,102,153,184,96,84,4,5674,	0,1412,20063,24336,25080,25800,25982,26182,
      my$chrom=$count_line_parts[0];
      my$start=$count_line_parts[1];
      $start=$start-1;
      my$end=$count_line_parts[2];
      my$fullcoords="$chrom:$start-$end";
      my$g="9";
      # get quantification
      my$num_counts=$count_line_parts[4];# still right
      if($num_counts>1){
        # get strand
        my$strand= $count_line_parts[3];
        # get annotation
        my$annot=$count_line_parts[8];
        # cleaning everything except first annotation
        if($annot=~/\,/){
          $annot=$`;
          #print "more than one annot!\n";
        }
        # use the next column if no annotation was found, else use next next
        if(!($annot=~/N/)){
          $annot=$count_line_parts[9];
          if(!($annot=~/N/)){
            $annot=$count_line_parts[10];# sometimes this is also the gene,so we are fine with that
          }
        }

        if((grep(/$fullcoords/,@allmapco))){
          $g=$mapping_hash{$fullcoords};
          #  print "found a match for $fullcoords\n";
        }
        my$line="$fullcoords\t$strand\t$samplename\t$num_counts\t$g\t$g\t$annot";
        $line=~s/\n//g;
        if(!($fullcoords=~/chrUn_gl/)){ # check for unsure coordinates
             print OU "$line\n";
        }
      }
  }
}
