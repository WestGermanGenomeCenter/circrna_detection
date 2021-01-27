#!/usr/bin/perl -w
use strict;
use Parallel::ForkManager;
use Getopt::Long qw(GetOptions);
# get the candidatelist_auto_all_sites.bed.csv file created with steptwo.pl
############################### example run
# perl ../../../auto_circs/auto_find_circ/matrixmaker-V4.pl --infile allsites_bedgroup_WNT.csv -outfile test6_mm4_matrixone_test1.tsv --logfile log.log --circ_bed_file ../../../auto_circs/auto_find_circ/circbase_known_circs.bed --gene_mapping_file ../../../auto_circs/auto_find_circ/genes_to_refseqID_nc_and_nr.tsv
##########################################
# made by Daniel R.

# input parameters: if not given, use default files:
####################################################
my $infile="in.tsv";
my$out_file="out.mat1";
my$logfile="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/logfile_auto.log";
my$circ_bed_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/circbase_known_circs.bed";
my$gene_mapping_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/genes_to_refseqID_nc_and_nr.tsv";
my$threads=12;
my$check_refseq=1;

GetOptions('circ_bed_file|c=s' => \$circ_bed_file,'gene_mapping_file|g=s' => \$gene_mapping_file,'logfile|l=s' => \$logfile,'threads|t=i' => \$threads,'infile|i=s' => \$infile,'outfile|o=s' => \$out_file,'check_refseq|chseqi=i'=>\$check_refseq) or warn "Using default parameters: $0 --from NAME\n";
chomp ($out_file,$infile,$logfile,$circ_bed_file,$gene_mapping_file,$threads,$check_refseq);# for files with gene names instead of NM** as annotation

# start with logfile
open(ER,'>>',$logfile)||die "$!";		# global logfile
my$start = time;
print ER "started MM_v4 with parameters:\ninfile:$infile\noutfile:$out_file\nchecking refseqids: $check_refseq\n";
print ER "circs_coords_mapping.bed:$circ_bed_file\ngene_mapping_file:$gene_mapping_file\nthreandN: $threads\nstarting...\n";
my$linfile= $infile;
chomp $linfile;
########################################################################### gene mapping file reading into hash %mapping
my%mapping=();
open(MA,$gene_mapping_file)|| die "$!";
my@allemappings= <MA>;
# each line now one array part
print ER "reading gene mapping...\n";
foreach my $mapline (@allemappings){
	chomp $mapline;
	#my $pid = $pm->start and next;
	if(!($mapline=~/^$/)){
	#	print "$mapline\n";
		my@slit=split(/\s+/,$mapline);
		my$genene=$slit[0];
		$genene =~ s/\s+//g; # remove emptieness
		my$nnum="";
		if(scalar(@slit)>1){
			$nnum=$slit[1];# will be key
			if($nnum=~/[NX]/){# we only want refseqids here, NR* or NM/NR/XR/XM
				$nnum =~ s/\s+//g;
				$mapping{"$nnum"}="$genene";
				#	print "mapping now $nnum to gene $genene\n";# refseqid to gene name
			}
		}
		$nnum="";
	}
}
close MA;

# candidatelist_auto_all_sites.bed.csv file created with steptwo.pl
print ER "reading input file $linfile ...\n";
# output file second argument adding coordinates
open(IN,$linfile)|| die "$!";
############################################ get samlenames into array, get coordinates and basic info into arrays
my@allelines= <IN>; #input file
my%sample_hash;
my%basic_info_hash;
my%alle_coords_hash;
my%alle_new_info_try_hash; # coords and strand to avoid redundant coords + strands with only different refseqid
DATA_IN:
for (my$i=0;$i<scalar(@allelines);$i++){
	my$line_o_o=$allelines[$i];	# current line
	if (!($line_o_o=~/coordinates/)){# ignore header if there , but should not be present at this stage anyway
		get_names($line_o_o);
		sub get_names{
			my$line= shift(@_);
			my@parts=split(/\t+/,$line);
			my$cord=$parts[0];
			my$strand=$parts[1];
			my$Refseqid=$parts[6];
			my$namesmale=$parts[2];
			if(!(exists($sample_hash{$namesmale}))){
				$sample_hash{$namesmale}=$i;
			}
			if(!(exists($alle_new_info_try_hash{"$cord"}))){ # check for coords redundancy
				$alle_new_info_try_hash{"$cord"}=$i;
				$basic_info_hash{"$Refseqid\t$i"}=$i;
				$alle_coords_hash{"$cord\t$strand"}=$i;
			}
		}
	}
}
# without tie to keep the hashes in order we need to sort the hash keys by their value- i.e the line number of the infile
# sort the hash keys based on their values (line number in infile )
my@sorted_by_val_allenames = sort{$sample_hash{$a} <=> $sample_hash{$b}} keys %sample_hash;# sort the keys according to value
my@allenames=@sorted_by_val_allenames;		# sample names
my@sorted_by_val_allecooords = sort{$alle_coords_hash{$a} <=> $alle_coords_hash{$b}} keys %alle_coords_hash;
my@allecooords=@sorted_by_val_allecooords; # coordinates
my@sorted_by_val_allebasicinfo = sort{$basic_info_hash{$a} <=> $basic_info_hash{$b}} keys %basic_info_hash;
my@allebasicinfo=@sorted_by_val_allebasicinfo;
# reading known circs file now since we have all  infile info we need for now
my%allinfoonesamplehash;
my$sampleout;
my%known_circs=();
open(CI,$circ_bed_file)|| die "$!";
my@alleci= <CI>;
################################################ gene mapping file reading into hashknown_circs
print ER "reading known circs...\n";
foreach my $circline (@alleci){			# fill a hash that is used later
	chomp $circline;
	if($circline=~/[a-z]/){
		my@slit=split(/\s+/,$circline);
		my$chr=$slit[0];
		my$cordst=$slit[1];
		my$cordnd=$slit[2];
		my$circname=$slit[3];
		my$fullcordmap="$chr:$cordst-$cordnd";
		chomp $fullcordmap;
		$known_circs{"$fullcordmap"}="$circname";
	}
}
close CI;
################### get all information from one sample into a hash , key is samplename and value is all information in one var
foreach my $samplenames (@allenames){
	$sampleout= `grep -w $samplenames $linfile`;	#
	$allinfoonesamplehash{"$samplenames"} = "$sampleout";
}
my@samples= keys %allinfoonesamplehash;
#print"should be all samplenames @samples\n";
my$outfile=$out_file;
chomp $outfile;
print ER "starting $threads threads, writing $outfile\n";
open(OU,">",$outfile)|| die "$!";
############################################# get stable header, build resizeable header for samples
print OU "coordinates\tstrand\tRefseqID\tGene\tknown_circ\tnum_samples_present\ttotal_sum_unique_counts\tqualities\tpresent_in_sample\t";
foreach my $sampls  (@allenames) {
	print OU "sample\t-unique_count\t-qualA\t-qualB\t"; # $sampls not in same order as below, need to change it
}
print OU "\n";
############################################# look for each circ in each sample and build a matrix
					# number of cores
my $pf = Parallel::ForkManager->new($threads);
my$ni=0;
our$count=0;
findc(\@allecooords);
sub findc{
	my@c=@{$_[0]};
  DATA_OUT:
  foreach my $circs (@allecooords){
  	$count ++;
		$pf->start and next DATA_OUT;
  	find_circ($circs);
  	sub find_circ {
    	my$circcand= shift(@_);
			$circcand=~s/\t[+-]{1}//;
			my$str=$&;# strand of circ candidate
			my$basicinfo=$allebasicinfo[$count -1]; # start at zero, while the counter is already at 1 in first iteration
			$basicinfo=~s/\t[0-9]{1,30}//;# remove the unique marker of maybe double circ information
			if($basicinfo=~/[A-z]/g){
				chomp $circcand;
				my$circn="unknown";
				chomp $basicinfo;
				my$presencething=""; # for each circ cand, add names of sapmles where it is present
				my$totalcounts=0;	# for each circ cand, add unique counts
				my$allquas=""; 		# for each sample, summarize qualities
				my$line="$circcand\t$basicinfo\t";
				$line=~s/\n//g;
				my$tolookup=$basicinfo;
				chomp $tolookup;
				my$allsamplelines="";
				my$allsamplehit=0;
				my$gene_name="";
				######### mapping circ to refseq, gene ###############
				if(exists($mapping{$tolookup})){
					my$geneo=$mapping{$tolookup};
					$line="$line\t$geneo";
					$gene_name=$geneo;
					#		print "found gene $gene_name for circ $tolookup in gene mapping hash\n";
				}
				else {
					$line="$line\tunkn";
					$gene_name="unkn";
				}
				if(exists($known_circs{$circcand})){
					$circn=$known_circs{$circcand};
				}
				else{
					$circn="unknown";
				}
				foreach my $single_sample (@allenames) {# looking in each sample for each circ
	  				my$allonesample= $allinfoonesamplehash{$single_sample};
        		if($allonesample=~/$circcand\t*.*\n/gi){### is the circ is found in sample### make this more explicit: rule out coords+1 at the end by \t??
							my$line_of_i=$&;
          		my$lineonesample=$line_of_i; #declare the interesting line
							$lineonesample=~s/\n//g;
							my@hit_line_parts=split(/\s+/,$lineonesample);
							# this line should be like chr1:117402185-117420649	+	trimmed.Lena_12__2	38	1	1 NR_00475
							my$hit_qunat=$hit_line_parts[3];
							my$strand=$hit_line_parts[1]; # should be the same as in $str # check
							my$hit_qualA=$hit_line_parts[4];
							my$hit_qualB=$hit_line_parts[5];
							my$ref_seq=$hit_line_parts[6];

							my$single_sample_hit_string="$single_sample\t$hit_qunat\t$hit_qualA\t$hit_qualB\t";
  	      		$allsamplelines="$allsamplelines$single_sample_hit_string";# attaching current result to all others of the same coords
							if($allsamplelines=~/chr/){
								warn "error in file: $allsamplelines should not include coordinates\nfull line: $line_of_i\n";
							}
							if($check_refseq){# for files with gene names instead of NM/NR/XR/XM
								if(($hit_qunat=~/[A-z]/)||($hit_qualA=~/[A-z]/)||($hit_qualB=~/[A-z]/)||(!($ref_seq=~/[NX][RMC]_[0-9]{3,11}/))){# sanity checks on the current line of interest
									warn "found mistake in line, no Refseqid?: $lineonesample\n";
								}
							}
	      			$presencething="$presencething-$single_sample";
	      			$allquas = "$allquas;$hit_qualA,$hit_qualB";
	      			$allquas =~s/\s+//g;
	      			$totalcounts=$totalcounts + $hit_qunat;
	      			$ni=$totalcounts;
	      			$allsamplehit++;
	  				}
        		else{# new: if circ not seen in current sample
							chomp $single_sample;
	      			$allsamplelines="$allsamplelines$single_sample\t0\t0\t0\t";
	  				}
  			}
				$basicinfo=~s/\n//g;
				$gene_name=~s/\n//g;
				if(($circcand=~/\:/)&&($presencething=~/[A-z]/)&&($circcand=~/^chr/)&&($str=~/[+-]/)&&(!($allquas=~/[A-z]/gi))){
  					my$linestring="$circcand\t$str\t$basicinfo\t$gene_name\t$circn\t$allsamplehit\t$ni\t$allquas\t$presencething\t$allsamplelines\n";
	  				$linestring  =~s/\t\t/\t/g;
	  				print OU $linestring;
	  				$linestring="";
						$gene_name="";
				}
				else{			# in case something with the line is wrong
	  				warn "error in line: circand is $circcand \n basicinfo is $basicinfo \n and presencething is $presencething\n";
				}
  		}
  	}
  	$pf->finish;
  }
}# findc end
$pf->wait_all_children;
my$end=time;
my$used_mins=($end-$start)/60;
print ER "done with matrix creation of file $outfile with input $linfile\nBuilding the matrix took $used_mins minutes\n";


###### descriptions ############################################################################################################
################################################################################################################################
#	matrixmaker.pl for find_circ (versions 1-4)
#		- takes the outfiles from steptwo (processed.csv/tsv)
#		- for more than one sample (more useful matrix) cat the .processed.csv files into one big file
#		- adds a little relevant information to each candidate
#		- needs an output filename- it will output in .tsv format to be readable for matrixtwo.pl (look at the bottom of this file
#			 for example lines for each relevant file)
#		- tracks time of usage
#		- dumps errors into  global logfile
##########################################
#
############################################# starting- getting input vars


# infile:
#coordinates\tstrand\tsample\treads\tqualA\tqualB\trefseqID
#chr10:100260217-100262063	-	ak_055_test19	2	40	5	NM_001303404

# $outfile will be the matrix1 file like : coordinates	strand	RefseqID	Gene	known_circ	num_samples_present	total_sum_unique_counts	qualities	present_in_sample
#chr10:102683731-102685776	+	NM_001136123	SLF2	hsa_circ_0006654	17	90	,40;5,40;40,	-NPTH61_S1_L_2-NPTH62_S2_L_2-NPTH63_S3_L_2-

# $logfile will be where some info will be dropped for debugging and the current stage.

# $circ_bed file will be a bed file with known circs from circbase(circbase.org, downloads-> all->.bed )

# $gene_mapping_file is a file containing refseqids-> gene name mapping like : 	VN1R48P
#														CHST14	NM_130468
