#  For complete documentation for circs_snake please go to doc/ in this repo
circs_snake: a snakemake workflow for circular RNA detection and quantification
## how to get started

- get the code: git clone https://github.com/daaaaande/circs_snake
- personalize: nano config.yaml
- prepare: nano samples.tsv
- also create dirs needed: dcc_2_out, f_c2_out, cx2_out
- run: snakemake -n (remove the -n once everything looks OK)
