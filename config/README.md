# users manual to circs_snake
circs_snake is a multi-pipeline circRNA detection workflow from RNASeq data.

This readme is meant to help you, the user, to understand what circs_snake tries to do such that you can use /change this to your liking / environment.
For an first rough overview, lets look at a DAG of this pipeline with two input samples.

![Alt text](dag_before_vote.svg?raw=true "circs_snake DAG")

Here you can see that (starting from the top) we have four major "starting points":
1. the parental pipeline flow (starting with rule r01): does the vote, normalization and preparation steps
2. find_circ (starting with fc_b, fc_a is a rule unpacking .fastq.gz files if this is the given format)
3. DCC (starting with dcc_b, dcc_a is a rule unpacking .fastq.gz files if this is the given format)
4. CIRCexplorer1 (starting with cx_b, cx_a is a rule unpacking .fastq.gz files if this is the given format)

each of the pipelines in run twice here, since we have two input samples in this example. The exception is the parental pipeline, this part will be only run once for each dataset.
Another visualization of the same flow is below, making this a little more clear:

![Alt text](crude_flowchart_circs.svg?raw=true "circs_snake flowchart")

Here you can see what happens with the data:
First all three pipelines (find_circ, DCC, CIRCexplorer1) are run on each sample, resulting in one file for each sample for each pipeline. An example output file at this stage looks like this:

![Alt text](example_tsv_file.png?raw=true "head of example tsv out file")

These files are summarized in step r06a,b,c that result in a .mat1 file for each pipeline.
The columns in this fle are: circRNA coordinates, strand, samplename, detected quantity, quality, quality, refseq annotation
Annotation is added, data is summarized and results in a .mat2 file (r07a,b,c).
These pipeline-specific matrix2 files are then voted (circRNA coordinates are overlapped and filtered based on only 3/3 overlaps) and finally then normalized, resulting in three normalized and voted circRNA datafiles as the main output of this pipeline. An example output file is given with example_output_norm_voted_dcc_hg19.csv


## before you can run this
Before you will be able to run this workflow, you need to have:
- snakemake installed
- have the find_circ scripts from the officical website (http://circbase.org/cgi-bin/downloads.cgi, Custom scripts for finding circRNAs; unpack, edit find_circ_conf.yaml accordingly)
- installed DCC and CIRCexplorer1 (install or download, edit the config.yaml files accordingly)
- reference genome index built for STAR and Bowtie2, aswell as the reference genome in .fa and .gtf format (other annotation data is in the data/ dir, edit the config.yaml files accordingly)
- all other software dependencies should be handled by snakemake, see the env.yaml files
- the config.yaml files are for **my** specific deployment, yours **should** vary. Here you only need to change directories for each of the needed files / folders + you can change pipeline-specific parameters to your liking aswell. I attached hg19 and hg38 example config.yaml files to ease your adaption.

and thats it! an example of how to execute the pipeline is given in howtostart.sh, a cluster config example is given in cluster_config.yaml and an example samplesheet is given aswell (samples.tsv)


## the samplesheet and expected files

Given this as samples.tsv:
``` head samples.tsv
samples
"SRR3184300"
"SRR3184285"
```
the workflow expects:

SRR3184300_1.fastq and SRR3184300_2.fastq +
SRR3184385_1.fastq and SRR3184385_2.fastq
 in the root directory of this workflow: path/to/circs_snake/. <- put the .fastq files here
 the lane identifier is changeable in the config.yaml:
 ```
 lane_ident1:
  "_1"
 lane_ident2:
  "_2"
```

The workflow itself does create the needed .tsv file given two input fastq files in its root directory. you can also self-create this, see scripts/snake_infile_creator.pl. (parental Snakefile, rule r03 is where this would be created from a previously created .fastq file list, rule r02)



## how to start a typical circs_snake run:
- copy/past/move paired end, trimmed and QC'ed .fastq files into circ_snake/.
- check if the lane idetifier is correct in all config.yaml files (change this if needed)
- `snakemake` (for more options here see howtorun.sh)



## further reading

For documentation on each single step, please refer to the original pipeline documentation: https://gitlab.com/daaaaande/circs/-/blob/master/README.md
