# execute snake_circs- cold start
cir
qsub -q default -I exec/interacive_job_snake.sh
screen -dRR snake_testrun1
module load Miniconda/3_snakemake
module load Snakemake


~/miniconda3/bin/snakemake -n
~/miniconda3/bin/snakemake -p  --cluster-config cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l walltime={cluster.time} -l select={cluster.nodes}:ncpus={cluster.ncpus}:mem={cluster.mem}:arch={cluster.arch}" -j 100 --latency-wait 900 --use-conda --max-status-checks-per-second 1 --keep-going -n
~/miniconda3/bin/snakemake -p  --cluster-config cluster.json --cluster "qsub -A {cluster.account} -q {cluster.queue} -l walltime={cluster.time} -l select={cluster.nodes}:ncpus={cluster.ncpus}:mem={cluster.mem}:arch={cluster.arch}" -j 100 --latency-wait 900 --use-conda --max-status-checks-per-second 1 --keep-going









##################################################
################ old stuff########################
conda avitvate base
env| grep python
python


export PYTHONPATH=/software/python/3.6.5/ivybridge/
export PYTHONHOME=/software/python/3.6.5/ivybridge/

export PYTHONHOME=~/miniconda3/bin/python3

export CONDA_EXE=~/miniconda3/bin/conda
export CONDA_PREFIX_1=/home/daric102/miniconda3/bin/conda
export LIBPATH=/home/daric102/miniconda3//lib
export CONDA_PYTHON_EXE=/home/daric102/miniconda3/bin/python
export PATH=/home/daric102/miniconda3/bin:/home/daric102/miniconda3/condabin:/usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/opt/pbs/bin:/home/daric102/bin
export LD_LIBRARY_PATH=/home/daric102/miniconda3/lib/
export MANPATH=/home/daric102/miniconda3/share/man/

