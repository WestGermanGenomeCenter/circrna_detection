# snakefile for find_circ hg38
# needs its own configfile, use this for now
configfile: "f_c/find_circ_snake_config.yaml"
samplesfile = "samples.tsv"
#
####  trying to load the samplesfile into a list
import pandas as pd

samples_df = pd.read_table(samplesfile).set_index("samples", drop=False)
sample_names = list(samples_df['samples'])
print("will now execute find_circ hg38 on samples:")
print(sample_names)

#rule all:
#    input:
#        expand(config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names)





rule reformat_circs:
    input:
        annotated_circs=config['prefix'] + "/f_c2_out/run_{name}/{name}.circ_candidates_auto.bed.annotated"
    params:
        strict=config['strict'],
        perl_script_dir=config['perl_scripts']
    output:
        tsv_out=config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv"
    shell:
        "perl {params.perl_script_dir}/f_c_outreader.pl --i {input} --o {output} -strict {params.strict}"

#       my$err8 = `grep circ $sample_dir_out/run_$sample_name.sites.bed | grep -v chrM | python2.7 $scripts_dir/sum.py -2,3 | python2.7 $scripts_dir/scorethresh.py -16 1 | python2.7 $scripts_dir/scorethresh.py -15 2 | python2.7 $scripts_dir/scorethresh.py -14 2 | python2.7 $scripts_dir/scorethresh.py 7 2 | python2.7 $scripts_dir/scorethresh.py 8,9 35 | python2.7 $scripts_dir/scorethresh.py -17 100000 >$sample_dir_out/run_$sample_name.circ_candidates_auto.bed`;

rule annotate_candidates:
    input:
        filtered_circs=config['prefix'] + "/f_c2_out/run_{name}/{name}.circ_candidates_auto.bed"
    output:
        annotated_circs=config['prefix'] + "/f_c2_out/run_{name}/{name}.circ_candidates_auto.bed.annotated"
    params:
        bed_file=config['annotation_bed']
    shell:
        "bedtools window -a {input} -b {params.bed_file} -w 1 >{output}"

rule filter_fc:
    input:
        bed_ci=config['prefix'] + "/f_c2_out/run_{name}/{name}.sites.bed"
    params:
        sum_script=config['fc_scripts_dir'] + "/sum.py",
        score_script=config['fc_scripts_dir'] + "/scorethresh.py"
    output:
        filtered_circs= config['prefix'] + "/f_c2_out/run_{name}/{name}.circ_candidates_auto.bed"
    shell:
        "grep circ {input.bed_ci} | grep -v chrM | python2.7 {params.sum_script} -2,3 | python2.7 {params.score_script} -16 1 | python2.7 {params.score_script} -15 2 | python2.7 {params.score_script} -14 2 | python2.7 {params.score_script} 7 2 | python2.7 {params.score_script} 8,9 35 | python2.7 {params.score_script} -17 100000 >{output} "

# bowtie2 --reorder --mm --score-min=C,-15,0 -q -x $bt_ref -U $sample_dir_out/auto_anchors.qfa 2> $sample_dir_out/auto_bt2_secondpass.log | python2.7 $scripts_dir/find_circ.py -G $fastqs_files/ -p $sample_name -s $sample_dir_out/run_$sample_name.sites.log > $sample_dir_out/run_$sample_name.sites.bed 2> $sample_dir_out/run_$sample_name.sites.reads
rule align_anchors:
    input:
        anchors=config['prefix'] + "/f_c2_out/run_{name}/{name}.auto_anchors.qfa"
    output:
        bed_ci=config['prefix'] + "/f_c2_out/run_{name}/{name}.sites.bed",
        reads_ci=config['prefix'] + "/f_c2_out/run_{name}/{name}.sites.reads",
        log_bt=config['prefix'] + "/f_c2_out/run_{name}/{name}.secondpass.log",
        log_fc=config['prefix'] + "/f_c2_out/run_{name}/{name}.f_c_run_sites.log"
    params:
        fastq_files_dir=config['chrom_fastqs'],
        genome=config['bowtie2_ref'],
        sample_name="{name}",
        fc_main_script=config['fc_scripts_dir'] + "/find_circ.py"
    shell:
        "bowtie2 --reorder --mm --score-min=C,-15,0 -q -x {params.genome} -U {input} 2> {output.log_bt} | python2.7 {params.fc_main_script} -G {params.fastq_files_dir} -p {params.sample_name} -s {output.log_fc} > {output.bed_ci} 2> {output.reads_ci}"



rule get_anchors:
    input:
        unmapped=config['prefix'] + "/f_c2_out/run_{name}/{name}.unmapped.auto.bam"
    output:
        anchors=config['prefix'] + "/f_c2_out/run_{name}/{name}.auto_anchors.qfa"
    params:
        fc_script=config['fc_scripts_dir'] + "/unmapped2anchors.py"
    shell:
        "python 2.7 {params.fc_script} {input} >{output}"


rule prepare_remapping:
    input:
        sorted_bam=config['prefix'] + "/f_c2_out/run_{name}/{name}.auto.bam"
    output:
        unmapped=config['prefix'] + "/f_c2_out/run_{name}/{name}.unmapped.auto.bam"
    threads:
        config['threads']
    shell:
        "samtools view -@ {threads} -hf 4 {input} | samtools view -@ {threads} -Sb - > {output}"
#       my$err4 = `samtools view -@ $threads -hf 4 $sample_dir_out/$sample_name.auto.bam | samtools view -@ $threads -Sb - > $sample_dir_out/unmapped_auto.bam`;


rule sort_bam:
    input:
        bam=config['prefix'] + "/f_c2_out/run_{name}/temp.bam"
    output:
        sorted_bam=config['prefix'] + "/f_c2_out/run_{name}/{name}.auto.bam"
    threads:
        config['threads']
    shell:
        "samtools sort @ {threads} -O bam -o {output} {input}"
# samtools sort -@ $threads -O bam -o $sample_dir_out/$sample_name.auto.bam $sample_dir_out/temp.bam
# now all the samtools fun
rule create_bam:
    input:
        sam=config['prefix'] + "/f_c2_out/run_{name}/temp.sam"
    output:
        bam=config['prefix'] + "/f_c2_out/run_{name}/temp.bam"
    threads:
        config['threads']
    shell:
        "samtools view -@ {threads} -hbuS -o {output} {input}"

rule map_primary_bt2:
    input:
        dir_check=config['prefix'] + "/f_c2_out/run_{name}/chk.tx",
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq"
    params:
        genome=config['bowtie2_ref'],
        logfile=config['prefix'] + "/f_c2_out/run_{name}/bt2_first_align.log"
    output:
        sam=config['prefix'] + "/f_c2_out/run_{name}/temp.sam"
        # bowtie2 -p $threads --very-sensitive --mm --score-min=C,-15,0 -x $bt_ref -1 $infile_dir/$lineonefile -2 $infile_dir/$linetwofile 2> $sample_dir_out/firstpass.log >$sample_dir_out/temp.sam\n"
    threads:
        config['threads']
    shell:
        "bowtie2 -p {threads} --very-sensitive --mm --score-min=C,-15,0 -x {params.genome} -1 {input.lane_1} -2 {input.lane_2} 2>{params.logfile} >{output.sam}"



rule prepare_files_fc:
    # unzip if needed,
    # create output dir
    #eed_unpack=0
    # create config["samples"] between-files dir
    # go to processing dir
    input:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq",
    params:
        processing_dir= config['prefix'] + "/f_c2_out/run_{name}/"
    output:
        config['prefix'] + "/f_c2_out/run_{name}/chk.tx"

    shell:
        "mkdir -p {params.processing_dir} && cd {params.processing_dir} && touch {output}"



rule unzip_gz_fc:
    input:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq",
        #out_dir=config['prefix'] + "/"
    output:
        lane_1=config['prefix'] + "/{name}" + config["lane_ident1"] + ".fastq.gz",
        lane_2=config['prefix'] + "/{name}" + config["lane_ident2"] + ".fastq.gz",
    shell:
        "gunzip {input} "