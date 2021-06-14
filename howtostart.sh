# execute snake_circs- cold start
cir
qsub -q default -I exec/interacive_job_snake.sh
screen -dRR snake
module load Miniconda/3_snakemake
module load Snakemake/5.10.0


~/miniconda3/bin/snakemake -n
~/miniconda3/bin/snakemake -p  --cluster-config cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l walltime={cluster.time} -l select={cluster.nodes}:ncpus={cluster.ncpus}:mem={cluster.mem}:arch={cluster.arch}" -j 100 --latency-wait 90000 --use-conda --max-status-checks-per-second 1 --keep-going -n
~/miniconda3/bin/snakemake -p  --cluster-config cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l walltime={cluster.time} -l select={cluster.nodes}:ncpus={cluster.ncpus}:mem={cluster.mem}:arch={cluster.arch}" -j 100 --latency-wait 90000 --use-conda --max-status-checks-per-second 1 --keep-going
