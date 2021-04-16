# circs_snake
re-write of circs, for now only hg38 paired-end.

### how to get started
- get the code:  `git clone https://github.com/daaaaande/circs_snake`
- personalize: `nano config.yaml`
- prepare: `nano samples.tsv`
  - also create dirs needed: dcc_2_out, f_c2_out, cx2_out
- run: `snakemake -n` (remove the -n once everything looks OK)
