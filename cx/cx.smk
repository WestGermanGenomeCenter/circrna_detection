# cx2 exec smk file
configfile: "cx/cx_conf.yaml"
samplesfile = "samples.tsv"

####  trying to load the samplesfile into a list
import pandas as pd

samples_df = pd.read_table(samplesfile).set_index("samples", drop=False)
sample_names = list(samples_df['samples'])
print("will now execute cx hg38 on samples:")
print(sample_names)


#rule all: # extend this rule to include outfiles from the last step respectively
#  input:
#    expand(config['prefix'] + "/cx_2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names)


# perl $infile_dir/circexplorer1_auto/circexplorer1_out_reader.pl $exec_dir/cx2_out_$sample_name.tsv $exec_dir/parsed_cx2_out_$sample_name.tsv $sample_name
rule cleanup_cx2anot:
    input:
        cx2_out=config['prefix'] + "/cx_2_out/run_{name}" + "/cx2_out.tsv"
    params:
        sample_name="{name}",
        pl_exec=config['pl_script_dir'] + "/circexplorer1_out_reader.pl"
    output:
        processed_fin=config['prefix'] + "/cx_2_out/run_{name}" + "/processed_run_{name}.tsv"
    shell:
        "perl {params.pl_exec} {input} {output} {params.sample_name}"

# $cx2_command_exec -r $refseq_file -g $fastqs_files -b $sample_dir_out/back_spliced_junction.bed -o $sample_dir_out/$sample_name.processed.tsv\n
#  -r $bowtie_index_dir/hg_38_ref_fetched.txt -g $bowtie_index_dir/hg38_full.fa -b $exec_dir/back_spliced_junction.bed -o $exec_dir/cx2_out_$sample_name.tsv`;
rule annotate_circs_cx2:
    input:
        junc=config['prefix'] + "/cx_2_out/run_{name}" + "/back_spliced_junction.bed"
    params:
        refseq_file=config['refseq_anot_file'],
        annotate_cmd=config['cx_annotate_command'],
        fasta_files=config['fasta_ref']
    output:
        cx2_out=config['prefix'] + "/cx_2_out/run_{name}" + "/cx2_out.tsv"
    shell:
        "{params.annotate_cmd} -r {params.refseq_file} -g {params.fasta_files} -b {input} -o {output}"



# following the cx2 flow
rule run_cx2:
    input:
        chim_out=config['prefix'] + "/cx_2_out/run_{name}" + "/Chimeric.out.junction"
    params:
        cx_parse_command=config['cx_parse_command'],
        log=config['prefix'] + "/cx_2_out/run_{name}" + "/CIRCexplorer2_parse.log"
    output:
        junc=config['prefix'] + "/cx_2_out/run_{name}" + "/back_spliced_junction.bed"
    shell:
        "{params.cx_parse_command} {input.chim_out} >{params.log}"



rule align_primary_STAR_CX:
  input:
    lane_1= config['prefix'] + "/" + "{name}" + config["lane_ident1"] + ".fastq",
    lane_2= config['prefix'] + "/" + "{name}" + config["lane_ident2"] + ".fastq",
    genome= config["genome_STAR_hg38"],
    dir_check= config['prefix'] + "/cx_2_out/run_{name}/chk.tx"
  params:
    processing_dir= config['prefix'] + "/cx_2_out/run_{name}/",
    outfile_name="{name}"
  threads: 12

  output:
    config['prefix'] + "/cx_2_out/run_{name}" + "/Chimeric.out.junction"

  conda:
        "/gpfs/project/daric102/circs_hilbert_scratchgs/snakemake_tests/envs/hg38_STAR.yaml"
  shell:
        "STAR --runThreadN {threads} --genomeDir {input.genome} --readFilesIn {input.lane_1} {input.lane_2} --chimSegmentMin 10 --limitBAMsortRAM 512000000000 "
#

rule prepare_files_cx2:
    # unzip if needed,
    # create output dir
    #eed_unpack=0
    # create config["samples"] between-files dir
    # go to processing dir
    input:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq",
    params:
        processing_dir= config['prefix'] + "/cx_2_out/run_{name}/",
    output:
        config['prefix'] + "/cx_2_out/run_{name}/chk.tx"

    shell:
        "mkdir -p {params} && cd {params} && touch {output}"



rule unzip_gz_cx2:
    input:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq",
        #out_dir=config['prefix'] + "/"
    output:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq.gz",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq.gz",
    shell:
        "gunzip {input} "
