configfile: "config.yaml"
samplesfile = "samples.tsv"

import pandas as pd

samples_df = pd.read_table(samplesfile).set_index("samples", drop=False)
sample_names = list(samples_df['samples'])
print("will now execute DCC, f_c and CX2 hg38 on samples:")
print(sample_names)

rule all: # extend this rule to include outfiles from the last step respectively
  input:
    expand(config['prefix'] + "/cx2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
    expand(config['prefix'] + "/f_c2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names),
    expand(config['prefix'] + "/dcc_2_out/run_{name}" + "/processed_run_{name}.tsv",name=sample_names)






include:"cx/Snakefile"
include:"dcc/Snakefile"
include:"f_c/Snakefile"
