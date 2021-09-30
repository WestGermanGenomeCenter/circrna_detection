#/usr/bin/perl -w
use strict;
use List::MoreUtils qw(uniq); # used to get to a list of unique samplenames later
use Getopt::Long qw(GetOptions);
############################################################ usage
# perl ../../../auto_circs/auto_find_circ/matrixtwo_V4.pl --mirRNA_file ../../../auto_circs/auto_find_circ/miRNA_circRNA_ineractions.txt --circbank_coding_file ../../../auto_circs/auto_find_circ/circRNA_protein_coding_potential.txt --hallmark_mapping_file ../../../auto_circs/auto_find_circ/hallmark_genes.tsv --ensembl_file ../../../auto_circs/auto_find_circ/mart_export_ensembl_gene_desc.txt --logfile ../matrixmaker_4/log.log --mapping_script ../../../auto_circs/auto_find_circ/read_mapping.pl --infile ../matrixmaker_4/test6_mm4_matrixone_test1.tsv -outfile test1_matrixtwo_getopt.mat2
# perl ../../auto_circs/auto_find_circ/matrixtwo_V4.pl --h ~/work_enclave/auto_circs/auto_find_circ/hallmark_genes.tsv --e ../../auto_circs/auto_find_circ/mart_export_ensembl_gene_desc.txt --n ../../auto_circs/auto_find_circ/read_mapping.pl --i mat1_examplefile_hg19.mat1 --o out_mat2_try2_nocircbank_allinfo.mat2 --l log.log --m ../../auto_circs/auto_find_circ/miRNA_circRNA_ineractions.txt --c ../../auto_circs/auto_find_circ/circRNA_protein_coding_potential.txt --excl_cb 1
############################################################

# made by Daniel R

# infiles. these are the default files- you can add your own files via extra parameters or change the file path
my$circbank_mirRNA_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/pipelines/miRNA_circRNA_ineractions.txt";
my$circbank_coding_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/pipelines/circRNA_protein_coding_potential.txt";
my$hallmark_mapping_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/hallmark_genes.tsv"; # unusual mapping file, not one gene per line
my$ensembl_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/mart_export_ensembl_gene_desc.txt";
my$logfile="logfile_auto.log";
my$mapping_script="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/read_mapping.pl";
# read custom files if provided
my$exclude_circbank=0;# set to 1 if you want to exclude circBank data (spliced length,mm_9 circ and prob_coding ),0 (default) includes this data then into the output
my$linfile="in.tsv";
my$outfile="out.tsv";



GetOptions('mirRNA_file|m=s' => \$circbank_mirRNA_file,'circbank_coding_file|c=s' => \$circbank_coding_file,'hallmark_mapping_file|h=s' => \$hallmark_mapping_file,'ensembl_file|e=s' => \$ensembl_file,'logfile|l=s' => \$logfile,'mapping_script|n=s' => \$mapping_script,'infile|i=s' => \$linfile,'outfile|o=s' => \$outfile,'exclude_circbank|excl_cb=i' => \$exclude_circbank) or warn "Using default parameters! \n";

chomp ($circbank_mirRNA_file,$circbank_coding_file,$hallmark_mapping_file,$ensembl_file,$logfile,$mapping_script,$linfile,$outfile);


############ all params now set
open(ER,'>>',$logfile)||die "$!";		# global logfile
print ER "reading input file $linfile ...\n";
# output file second argument adding coordinates
open(IN,$linfile)|| die "$!\n no infile ";
my@allelines= <IN>;
#my$outfile=$ARGV[1];  # outfile
#chomp $outfile;

##############################################################################
# start of all mapping
##############################################################################
require $mapping_script; # module reading mapping file for additional information- can be ignored when useless


open(ER,'>>',$logfile)||die "$!";		# global logfile
my $start = time;
open (OUT ,">",$outfile)|| die "$!";
open(MA,$hallmark_mapping_file) || die "$!";
# uses subroutine map_file from read_mapping.pl
my%mart_info=map_file($ensembl_file,1,2,"\t");
my@mart_infos= keys %mart_info;
my@allemappings= <MA>;
my%mapping_hash=();  # mapping hash: gene name is key, hallmark mechanism is value?
# each line now one array part
#print "reading gene mapping...\n";
######################### mapping hallmark file ###########
foreach my $mapline (@allemappings){
# fill a hash that is used later
  chomp $mapline;
  my@mappingline_parts=split(/\s+/,$mapline);
  my$hallmarktype_full=shift @mappingline_parts;# getting hallmark description properly with some regex cleaning by shifting first array index
  $hallmarktype_full=~s/HALLMARK//;
  $hallmarktype_full=~s/\s+//;
  $hallmarktype_full=~s/_//;
  my$address= shift @mappingline_parts;
  # rest of the line is all hallmark genes, need to be cleaned up before saving
  foreach my $hallmg (@mappingline_parts){
    $hallmg=~s/\s+//;
    if($hallmg=~/[A-Z]/){  # checking for empty lines and gene names
      $mapping_hash{"$hallmg"}="$hallmarktype_full";
    }
  }
}
my@allehallmarkg=keys %mapping_hash;
my@all_hm_genes= values %mapping_hash;
############################################################################
# arrays in use
my@sampleuniqc=();# positions of unique count columns
my@samplenames=();# names of all detected samples
my@uniqcounts=(); # where maybe all unique counts will be added into a two-dimensional array?
my@headers=();    # headers with relevant information for each circ candidates
my@alluniques=();# empty yet



# to be filled with circbank info
my%coords_to_circ_hash=();
my%coords_length_mouse_info_hash=();
my%circ_to_prob_hash=();

if(!($exclude_circbank)){# set parameter to 0 means we include the data, otherwise we ignore all circbank information
  open(MI,$circbank_mirRNA_file) || die "$!";
  my@allintermirrna=<MI>;



  #my%coords_mouse_cons=()
  ######################### mapping circbank file ###########

  foreach my $mir_cb_line (@allintermirrna){
  	my$mouse=0;
  	my$full_mouse_coords="Na";
  	my@parts_filem=split(/\s+/,$mir_cb_line);
  	# build into coordinates
  	my$chromosome=$parts_filem[2];
  	my$strt=$parts_filem[3];
  	my$nd=$parts_filem[4];
  	my$human_circrna_name=$parts_filem[0]; # we can map to this so we do not need coordinates- we match with circ name
  	my$strnd_map_circbank=$parts_filem[5]; # for later checking the strand
  	my$spliced_length=$parts_filem[6];# maybe use later
  	my$full_cordss="$chromosome:$strt-$nd";# this can be key or value now
  	my$circ_name_conv=$parts_filem[1];
  	$coords_to_circ_hash{"$circ_name_conv"}="$full_cordss";# circ name is key, coords are value
  	my$whole_line=join(/_/,@parts_filem);
  	# mouse circRNA detection
  	if($whole_line=~/mmu_circ/){
  		my$mouse_cords=$';
  		$full_mouse_coords="mmu_circ_$mouse_cords"
  	}
  	# key: coordinates , value= spliced_length and mouse coords if found
  	$coords_length_mouse_info_hash{$human_circrna_name}="$spliced_length\t$full_mouse_coords";
  }
  ######################### mapping circbank  coding potential file ###########

  open(CO,$circbank_coding_file) || die "$!";
  my@allincoding=<CO>;

  foreach my $codingcb_line (@allincoding){
  	my@parts_fileco=split(/\s+/,$codingcb_line);
  	my$circ_name=$parts_fileco[0];
  	my$coding_prob=$parts_fileco[5];
  	$circ_to_prob_hash{"$circ_name"}=$coding_prob;
  }
}
##############################################################################
# end of all mapping
##############################################################################
# default non-matching descriptions
my$marti="NaN";
my$hallm="none\t";
# starting to read inputfile line by line
for (my $var = 0; $var < scalar(@allelines); $var++) {
  my$longline=$allelines[$var];
  if ($var > 0) {   # header of a .mat1 file is ALWAYS ONLY THE FIRST LINE
    # getting relevant information for each circ candidate ...
    my@lineparts=split(/\t/,$longline);
    my$coords=$lineparts[0];
    my$strand=$lineparts[1];
    my$refseqID=$lineparts[2];
    my$gene=$lineparts[3];
    my$circn=$lineparts[4];

    # if refseq has NM/R in front, everything is good. if it has not, iths the gene name and we do not know the RefseqID

    if(!($refseqID=~/[NX][MRC]_[0-9]{3,11}/)){# now also allowing predicted mrna/ncrna
          $gene=$refseqID;
   }

    ############
    # starting to add information from all mapping files based on circRNA name, gene or coordinates
    # adding hallmark gene type
    ############
    $hallm="none\t";# default value if none found, will only be overwritten once a matching hallmark match is found
    if(grep(/$gene/,@allehallmarkg)){
    # find hallmark class  and add to matrix file
      if($mapping_hash{$gene}=~/[A-Z]/){
        $hallm=$mapping_hash{$gene};
      }
    }

    if(grep(/$gene?/,@mart_infos)){              # mart mapping
      # gene has information available to it
      $marti=$mart_info{$gene};
      $marti=~s/\[.*\]//g;
	$marti=~s/\ /_/g;
    }
    else{
	$marti="NaN";
    }
    # check for empty mart information
    if(!($marti=~/[A-z]/gi)){
      $marti="NaN";
    }
    # integration of circbank data
    if(!($exclude_circbank)){
      my$prob_coding="nA";
      if($circ_to_prob_hash{"$circn"}=~/[0-9]/){# get the coding probability based on hssa_circ_name in both files
  		  $prob_coding=$circ_to_prob_hash{"$circn"};
      }
      # mouse and length data
      my$length_and_mouse_circ="nA\tnone";
      if($coords_length_mouse_info_hash{"$circn"}=~/\W+/){
  		    $length_and_mouse_circ=$coords_length_mouse_info_hash{"$circn"};
      }
      push(@headers,"$coords\t$strand\t$refseqID\t$gene\t$circn\t$hallm\t$marti\t$length_and_mouse_circ\t$prob_coding");# header into array

    }
    else{# ecluding these columns
    push(@headers,"$coords\t$strand\t$refseqID\t$gene\t$circn\t$hallm\t$marti");# header into array

    }
    my$e=0;
    my$allthings="";
    foreach my$samplepos (@sampleuniqc){
      $e++; # second coordinate for two-dimensional array of all unique counts
      my$samplename = $lineparts[$samplepos-1];
      push(@uniqcounts,$lineparts[$samplepos]); #the unique
      #  print "found a sample= $samplename \nand its counts are $lineparts[$samplepos] circrna of interest is $lineparts[4]\n";
      push (@samplenames,$samplename);
      ## do the magic and find all unique counts for each sample for each circrna candindate, get all this into a string and then print all that out later
      $allthings="$allthings\t$lineparts[$samplepos]";
      #$alluniques[$var][$e]=$lineparts[$samplepos];# two-dimensional array
    }
    push(@alluniques,$allthings);# all unique count positions in .mat1 file
    ##############
    # output line finished at this point
    ##############
  }
  else{
  # header bar only- catch columns with samplenames in it and their position
    my@wholeheader=split(/\t/,$longline);
    my $i=0;
    foreach my $headername (@wholeheader){
      $i++;
      if ($headername=~/sample/) {
        if(!($headername=~/\_sample/)){
          push (@sampleuniqc,$i);  #get positions of sample ids
        }
      }
    }
}# header column number collection finished
}



# actual file creation:
# header line in output file ...
my@uniques= uniq @samplenames;
if(!($exclude_circbank)){
  print OUT"coordinates\tstrand\trefseqid\tgene\tcircn\thallm\tbiom_desc\tspliced_length\tmm9_circ\tprob_coding\t";
}
else{
  print OUT"coordinates\tstrand\trefseqid\tgene\tcircn\thallm\tbiom_desc\t";
}
foreach my $sampl (@uniques){
  print OUT"$sampl\t";
}
print OUT "\n";
# now the real content, cleaning it from junk and then printing it
for (my $v = 0; $v < scalar(@headers); $v++) {
  my$outline="$headers[$v]$alluniques[$v]";
  $outline=~s/\t\t+/\t/g;
  $outline=~s/\t\s+/\t/g;
  $outline=~s/\s+\t/\t/g;
  print OUT "$outline\n";
}

##### explanations #########################
############################################

# additional information from circBank: coding potential + mirRNA interaction sites

# mapping those files

# example line mirRNA file:
#id	circRNA_id	chr	start	end	strand	spliced_seq_length	annotation	best_transcript	gene_symbol	mm9_circRNA_id
#hsa_circA1CF_001	hsa_circ_0018410	chr10	52575765	52623840	-	1470	ANNOTATED, CDS, coding, INTERNAL, OVCODE, OVERLAPTX, OVEXON, UTR5	NM_001198819	A1CF	#N/A
#hsa_circA1CF_002	hsa_circ_0018409	chr10	52559168	52573822	-	7964	ANNOTATED, CDS, coding, OVCODE, OVERLAPTX, OVEXON, UTR3	NM_001198819	A1CF	#N/A
#hsa_circA2ML1_001	hsa_circ_0025378	chr12	8976315	8988935	+	482	ANNOTATED, CDS, coding, INTERNAL, OVCODE, OVERLAPTX, OVEXON	NM_144670	A2ML1	#N/A
#hsa_circA2ML1_002	hsa_circ_0025379	chr12	8990035	8998818	+	955	ANNOTATED, CDS, coding, INTERNAL, OVCODE, OVEXON	NM_144670	A2ML1	#N/A
#hsa_circA2ML1_003	hsa_circ_0025380	chr12	8993964	9001510	+	948	ANNOTATED, CDS, coding, INTERNAL, OVCODE, OVEXON	NM_144670	A2ML1	#N/A
# since there is no real interaction prediction, we here do the circbank coordinates to id conversion

# in the second file we can then get from id to coding probability

#			matrixtwo.pl
#			- needs an infile -> the correct infile format is made by matrixmaker.pl as ouput, can directly given to matrixtwo
#			- adds additional biologic information, but also removes information from the first matrix that is not used in the heatmap that will be created with the output from this script
#			- needs output file name , outputs in a .tsv file format
#     - will run in the dir where it was started
#			- will output errors into ER file : /home/daniel/logfile_auto.log, can be changed

############## example input :
#coordinates	strand	RefseqID	Gene	known_circ	num_samples_present	total_sum_unique_counts	qualities	present_in_sample	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB	sample	-unique_count	-qualA	-qualB
#chr10:102683731-102685776	+	NM_001136123	SLF2	hsa_circ_0006654	16	88	,6;40,40;40,40;5,40;5,6;40,6;40,6;40,40;40,40;5,40;40,40;40,40;40,40;40,40;5,40;40,40;40	-run_hal01_test1a-run_697_r_test1a-run_hal01_r_test1a-run_hal01_r_test1c-run_hal01_test1c-run_hal01_test1d-run_hal01_test1b-run_697_r_test1c-run_hal01_r_test1b-run_697_test1c-run_697_r_test1b-run_697_test1e-run_697_r_test1d-run_hal01_r_test1d-run_697_test1a-run_697_test1d	run_hal01_test1a	5	6	40	run_697_r_test1a	7	40	40	run_hal01_r_test1a	7	40	5	run_hal01_r_test1c	7	40	5	run_hal01_test1c	5	6	40	run_hal01_test1d	5	6	40	run_hal01_test1b	5	6	40	run_697_r_test1c	7	40	40	run_hal01_r_test1b	7	40	5	run_697_test1c	3	40	40	run_697_r_test1b	7	40	40	run_697_test1e	3	40	40	run_697_r_test1d	7	40	40	run_hal01_r_test1d	7	40	5	run_697_test1a	3	40	40	run_697_test1d	3	40	40
#chr10:102683734-102685776	+	NM_001136123	SLF2	unknown	4	8	,40;40,40;40,40;40,40;40	-run_697_r_test1a-run_697_r_test1c-run_697_r_test1b-run_697_r_test1d	run_hal01_test1a	0	0	0	run_697_r_test1a	2	40	40	run_hal01_r_test1a	0	0	0	run_hal01_r_test1c	0	0	0	run_hal01_test1c	0	0	0	run_hal01_test1d	0	0	0	run_hal01_test1b	0	0	0	run_697_r_test1c	2	40	40	run_hal01_r_test1b	0	0	0	run_697_test1c	0	0	0	run_697_r_test1b	2	40	40	run_697_test1e	0	0	0	run_697_r_test1d	2	40	40	run_hal01_r_test1d	0	0	0	run_697_test1a	0	0	0	run_697_test1d	0	0	0
#
############# example output :
#
#coordinates	refseqid	gene	circn	hallm	biom_desc	run_hal01_test1a	run_697_r_test1a	run_hal01_r_test1a	run_hal01_r_test1c	run_hal01_test1c	run_hal01_test1d	run_hal01_test1b	run_697_r_test1c	run_hal01_r_test1b	run_697_test1c	run_697_r_test1b	run_697_test1e	run_697_r_test1d	run_hal01_r_test1d	run_697_test1a	run_697_test1d
#chr10:102683731-102685776	NM_001136123	SLF2	hsa_circ_0006654	none	SMC5-SMC6_complex_localization_factor_2_	5	7	7	7	5	5	5	7	7	3	7	3	7	7	3	3
#chr10:102683734-102685776	NM_001136123	SLF2	unknown	none	SMC5-SMC6_complex_localization_factor_2_	0	2	0	0	0	0	0	2	0	0	2	0	2	0	0	0
#
#
#
#
########################################################################### gene mapping file reading into hash %mapping

##exapmle line : HALLMARK_TNFA_SIGNALING_VIA_NFKB	http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_TNFA_SIGNALING_VIA_NFKB	JUNB	CXCL2	ATF3	NFKBIA	TNFAIP3	PTGS2	CXCL1	IER3	CD83	CCL20	CXCL3	MAFF	NFKB2	TNFAIP2	HBEGF	KLF6	BIRC3	PLAUR	ZFP36	ICAM1	JUN	EGR3	IL1B	BCL2A1	PPP1R15A	ZC3H12A	SOD2	NR4A2	IL1A	RELB	TRAF1	BTG2	DUSP1	MAP3K8	ETS2	F3	SDC4	EGR1	IL6	TNF	KDM6B	NFKB1	LIF	PTX3	FOSL1	NR4A1	JAG1	CCL4	GCH1	CCL2	RCAN1	DUSP2	EHD1	IER2	REL	CFLAR	RIPK2	NFKBIE	NR4A3	PHLDA1	#IER5	TNFSF9	GEM	GADD45A	CXCL10	PLK2	BHLHE40	EGR2	SOCS3	SLC2A6	PTGER4	DUSP5	SERPINB2	NFIL3	SERPINE1	TRIB1	TIPARP	RELA	BIRC2	CXCL6	LITAF	TNFAIP6	CD44	INHBA	PLAU	MYC	TNFRSF9	SGK1	TNIP1	NAMPT	FOSL2	PNRC1	ID2	CD69	IL7R	EFNA1	PHLDA2	PFKFB3	CCL5	YRDC	IFNGR2	SQSTM1	BTG3	GADD45B	KYNU	G0S2	BTG1	MCL1	VEGFA	MAP2K3	CDKN1A	CYR61	TANK	IFIT2	IL18	TUBB2A	IRF1	FOS	OLR1	RHOB	AREG	NINJ1	ZBTB10	PPAP2B	KLF4	CXCL11	SAT1	CSF1	GPR183	PMEPA1	PTPRE	TLR2	CXCR7	KLF10	MARCKS	#LAMB3	CEBPB	TRIP10	F2RL1	KLF9	LDLR	TGIF1	RNF19B	DRAM1	B4GALT1	DNAJB4	CSF2	PDE4B
