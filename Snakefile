configfile: "config.yaml"
samplesfile = "samples.tsv"
runfile = "run_conf.yaml"
import pandas as pd

samples_df = pd.read_table(samplesfile).set_index("samples", drop=False)
sample_names = list(samples_df['samples'])
#print("will now execute DCC, f_c and CX2 hg38 on samples:")
#print(sample_names)






# mm2 creation :
#perl /gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/matrixtwo_V4.pl --i f_c2_out/test_array_mb30/fc2_matrix_infile_3.mat1 --o f_c2_out/test_array_mb30/fc2_matrix_infile_3.mat2 --excl_cb 1
#perl /gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/matrixtwo_V4.pl --i dcc_2_out/test_array_mb30/dcc2_matrix_infile_3.mat1 --o dcc_2_out/test_array_mb30/dcc2_matrix_infile_3.mat2 --excl_cb 1
#perl /gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/matrixtwo_V4.pl --i cx2_out/test_array_mb30/cx2_matrix_infile_3.mat1 --o cx2_out/test_array_mb30/cx2_matrix_infile_3.mat2 --excl_cb 1

# file paths to change:
#my$circbank_mirRNA_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/pipelines/miRNA_circRNA_ineractions.txt";
#my$circbank_coding_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/pipelines/circRNA_protein_coding_potential.txt";
#my$hallmark_mapping_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/hallmark_genes.tsv"; # unusual mapping file, not one gene per line
#my$ensembl_file="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/mart_export_ensembl_gene_desc.txt";
#my$logfile="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/logfile_auto.log";
#my$mapping_script="/gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/read_mapping.pl";

# parameter names
#GetOptions('mirRNA_file|m=s' => \$circbank_mirRNA_file,'circbank_coding_file|c=s' => \$circbank_coding_file,'hallmark_mapping_file|h=s' => \$hallmark_mapping_file,'ensembl_file|e=s' => \$ensembl_file,'logfile|l=s' => \$logfile,'mapping_script|n=s' => \$mapping_script,'infile|i=s' => \$linfile,'outfile|o=s' => \$outfile,'exclude_circbank|excl_cb=i' => \$exclude_circbank) or warn "Using default parameters! \n";


# from here on the post_run guide need to be implemented- try first to get this to run so far
# do 3 rules that do exactly the same , one job for each pipeline outfile


# to make :
# perl /gpfs/project/daric102/circs_hilbert_scratchgs/repo/circs/auto_find_circ/matrixtwo_V4.pl --i f_c2_out/test_array_mb30/fc2_matrix_infile_3.mat1 --o f_c2_out/test_array_mb30/fc2_matrix_infile_3.mat2 --excl_cb 1

rule all:
    input:
        cx2_mat2=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2.mat2",
        dc2_mat2=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2.mat2",
        fc2_mat2=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2.mat2"




rule _r07c_run_matrix2_cx:
    input:
        cx_hg38_mat1=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.mat1"
    params:
        perl_script_mm2=config['mm2_script'],
        micrornas_file=config['micrornas_file'],
        coding_circnas_file=config['circbank_coding_file'],
        hallmarks_file=config['hallmarks_file'],
        ensembl_file=config['ensembl_file'],
        mapping_script=config['mapping_script']
    output:
        cx2_mat2=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2.mat2"
    shell:
        "perl {params.perl_script_mm2} --i {input} --o {output} --m {params.micrornas_file} --c {params.coding_circnas_file} --h {params.hallmarks_file} --e {params.ensembl_file} --n {params.mapping_script} --e {params.ensembl_file} --excl_cb 1"



rule _r07b_run_matrix2_dc:
    input:
        dc_hg38_mat1=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.mat1"
    params:
        perl_script_mm2=config['mm2_script'],
        micrornas_file=config['micrornas_file'],
        coding_circnas_file=config['circbank_coding_file'],
        hallmarks_file=config['hallmarks_file'],
        ensembl_file=config['ensembl_file'],
        mapping_script=config['mapping_script']
    output:
        dc2_mat2=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2.mat2"
    shell:
        "perl {params.perl_script_mm2} --i {input} --o {output} --m {params.micrornas_file} --c {params.coding_circnas_file} --h {params.hallmarks_file} --e {params.ensembl_file} --n {params.mapping_script} --e {params.ensembl_file} --excl_cb 1"



rule _r07a_run_matrix2_fc:
    input:
        fc_hg38_mat1=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.mat1"
    params:
        perl_script_mm2=config['mm2_script'],
        micrornas_file=config['micrornas_file'],
        coding_circnas_file=config['circbank_coding_file'],
        hallmarks_file=config['hallmarks_file'],
        ensembl_file=config['ensembl_file'],
        mapping_script=config['mapping_script']
    output:
        fc2_mat2=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2.mat2"
    shell:
        "perl {params.perl_script_mm2} --i {input} --o {output} --m {params.micrornas_file} --c {params.coding_circnas_file} --h {params.hallmarks_file} --e {params.ensembl_file} --n {params.mapping_script} --e {params.ensembl_file} --excl_cb 1"


rule _r06c_run_matrixmaker_cx:
    input:
        all_cx2_out_catted=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.tx"
    params:
        perl_script_m1=config['mm1_script'],
        annotation_file_m1=config['mm1_refseq_file'],
        circs_bed_file=config['mm1_circ_bedfile']
    output:
        cx_hg38_mat1=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.mat1"
    shell:
        "perl {params.perl_script_m1} --i {input} --o {output} --c {params.circs_bed_file} --g {params.annotation_file_m1}"



rule _r06b_run_matrixmaker_dc:
    input:
        all_dcc2_out_catted=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.tx"
    params:
        perl_script_m1=config['mm1_script'],
        annotation_file_m1=config['mm1_refseq_file'],
        circs_bed_file=config['mm1_circ_bedfile']
    output:
        dc_hg38_mat1=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.mat1"
    shell:
        "perl {params.perl_script_m1} --i {input} --o {output} --c {params.circs_bed_file} --g {params.annotation_file_m1}"




rule _r06a_run_matrixmaker_fc:
    input:
        all_fc2_out_catted=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.tx"
    params:
        perl_script_m1=config['mm1_script'],
        annotation_file_m1=config['mm1_refseq_file'],
        circs_bed_file=config['mm1_circ_bedfile']
    output:
        fc_hg38_mat1=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.mat1"
    shell:
        "perl {params.perl_script_m1} --i {input} --o {output} --c {params.circs_bed_file} --g {params.annotation_file_m1}"





rule _r05c_collect_tsvs_cx: # extend this rule to include outfiles from the last step respectively
    input:
        cx2_outfiles=expand(config['prefix'] + "/cx2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        #fc2_outfiles=expand(config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        #dc2_outfiles=expand(config['prefix'] + "/dcc_2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        reads_per_samplefile=config["prefix"]+"/reads_per_sample_"+config['run_name']+".tsv" # to get the preparations to be executed
    params:
        run_name=config['run_name']
    output:
        #all_fc2_out_catted=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.tx",
        all_cx2_out_catted=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.tx"
        #all_dcc2_out_catted=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.tx"
    shell:
        "cat {input.cx2_outfiles} >{output.all_cx2_out_catted}"

rule _r05b_collect_tsvs_dc: # extend this rule to include outfiles from the last step respectively
    input:
        #cx2_outfiles=expand(config['prefix'] + "/cx2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        #fc2_outfiles=expand(config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        dc2_outfiles=expand(config['prefix'] + "/dcc_2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        reads_per_samplefile=config["prefix"]+"/reads_per_sample_"+config['run_name']+".tsv" # to get the preparations to be executed
    params:
        run_name=config['run_name']
    output:
    #    all_fc2_out_catted=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.tx",
    #    all_cx2_out_catted=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.tx",
        all_dcc2_out_catted=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.tx"
    shell:
        "cat {input.dc2_outfiles} >{output.all_dcc2_out_catted}"

rule _r05a_collect_tsvs_fc: # extend this rule to include outfiles from the last step respectively
    input:
        #cx2_outfiles=expand(config['prefix'] + "/cx2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        fc2_outfiles=expand(config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        #dc2_outfiles=expand(config['prefix'] + "/dcc_2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
        reads_per_samplefile=config["prefix"]+"/reads_per_sample_"+config['run_name']+".tsv" # to get the preparations to be executed
    params:
        run_name=config['run_name']
    output:
        all_fc2_out_catted=config['prefix']+"/f_c2_out/"+config['run_name']+"/all_"+config['run_name']+"_fc2_tsvs.tx",
        #all_cx2_out_catted=config['prefix']+"/cx2_out/"+config['run_name']+"/all_"+config['run_name']+"_cx2_tsvs.tx",
        #all_dcc2_out_catted=config['prefix']+"/dcc_2_out/"+config['run_name']+"/all_"+config['run_name']+"_dcc2_tsvs.tx"
    shell:
        "cat {input.fc2_outfiles} >{output.all_fc2_out_catted}"

# start with collecting reads information, utilize perl perl_scripts_di


# normalization of matrixtwo values

# matrixtwo execution x 3

# matrixmaker execution x3

# cleanup/correct .tsv files

# cat the tsvs into one

# collect all outfiles in dir for each pipeline






# execute 3x main pipeline





# make dir where all outfiles will be collected

# move all .fastqs to /repo/circs dir

# perl $perl_scripts_dir/fastqs_to_reads_per_sample.pl --i reads_infile_$proj_name.tx --m $paired_single_end_param >reads_per_sample_$proj_name.tsv
rule _r04_count_reads_fastqs:
    # make another new script
    input:
        infile=config["prefix"]+"/samples.tsv",
        fastq_list_file=config["prefix"]+"/fastq_infiles_list.tx"
    params:
        run_name=config['run_name'],
        perl_script2=config["perl_script_dir"]+ "/fastq_list_to_reads_per_sample.pl",
        lane_1ident=config["lane_ident1"]
    output:
        reads_per_samplefile=config["prefix"]+"/reads_per_sample_"+config['run_name']+".tsv"
    shell:
        "perl {params.perl_script2} --l1 {params.lane_1ident} --i fastq_list_file >{output}"


# perl $perl_scripts_dir/fastq_list_to_infile_circs_ext.pl --i fastq_list_$proj_name.tx --g $proj_name --l1 $lane_diff --l2 $lane_to_id --m $paired_single_end_param >infile_$proj_name.tx
rule _r03_create_infile:
    input:
        fastq_list_file=config["prefix"]+"/fastq_infiles_list.tx"
    params:
        perl_script=config["perl_script_dir"]+ "/snake_infile_creator.pl",
        run_name=config['run_name'],
        lane_1ident=config["lane_ident1"],
        lane_2ident=config["lane_ident2"]
    output:
        infile=config["prefix"]+"/samples.tsv" # workaround for now, delete the old infile before a new one get created
    shell:
        " rm {output} && perl {params.perl_script} --i {input} --l1 {params.lane_1ident} --l2 {params.lane_2ident} >{output}"
# ls -1f *.fastq>fastq_list_$proj_name.tx

# create dir where output files get copied over
rule _r02_create_fastq_list:
    input:
        fc_2_chk_to_create=config["prefix"]+"/f_c2_out/"+config['run_name']+"/chk.tx",
        cx_2_chk_to_create=config["prefix"]+"/cx2_out/"+config['run_name']+"/chk.tx",
        dc_2_chk_to_create=config["prefix"]+"/dcc_2_out/"+config['run_name']+"/chk.tx"
    params:
        dir_to_cd_to=config["prefix"]
    output:
        fastq_list_file=config["prefix"]+"/fastq_infiles_list.tx"
    shell:
        " cd {params.dir_to_cd_to} && ls -f1 *.fastq >{output}"


rule _r01_prepare_hg38_dirs:
    input:
        run_file="run_conf.yaml" # this might not work
    params:
        fc_2_dir_to_create=config["prefix"]+"/f_c2_out/"+config['run_name'],
        cx_2_dir_to_create=config["prefix"]+"/cx2_out/"+config['run_name'],
        dc_2_dir_to_create=config["prefix"]+"/dcc_2_out/"+config['run_name'],
    output:
        # 3x chk.tx files, as before
        fc_2_chk_to_create=config["prefix"]+"/f_c2_out/"+config['run_name']+"/chk.tx",
        cx_2_chk_to_create=config["prefix"]+"/cx2_out/"+config['run_name']+"/chk.tx",
        dc_2_chk_to_create=config["prefix"]+"/dcc_2_out/"+config['run_name']+"/chk.tx"
    shell:
        "mkdir -p {params} && touch {output}"


include:"cx/Snakefile"
include:"dcc/Snakefile"
include:"f_c/Snakefile"
