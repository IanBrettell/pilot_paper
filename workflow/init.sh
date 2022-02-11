# For snakemake

module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
cd /hps/software/users/birney/ian/repos/pilot_paper
conda activate snakemake_6.13.1
snakemake \
  --jobs 5000 \
  --latency-wait 100 \
  --cluster-config config/cluster.yaml \
  --cluster 'bsub -g /snakemake_bgenie -J {cluster.name} -q {cluster.queue} -n {cluster.n} -M {cluster.memory} -o {cluster.outfile}' \
  --keep-going \
  --rerun-incomplete \
  --use-conda \
  --use-singularity \
  -s workflow/Snakefile \
  -p

#######################
# Containers
#######################

# Open CV (python)
OPENCVCONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/opencv_4.5.1.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $OPENCVCONT \
    workflow/envs/opencv_4.5.1.def


# idtrackerai

bsub -M 20000 -q gpu -gpu "num=1:gmem=1000" -Is bash
module load cuda-10.0.130-gcc-9.3.0-336gifa
conda create -n idtrackerai python=3.6
conda activate idtrackerai
pip install idtrackerai[gpu]
# This connects the blobs with ~1.1 iterations/second
# This doesn't change by increasing `num` (from 1 to 2)
# Google Colab does it at ~10 iterations/second

CONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/idtrackerai.sif
singularity build --remote \
    $CONT \
    workflow/envs/20211125_idtrackerai.def

# For R

## Set container location
CONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/R_4.1.0.sif

## Build Rocker container
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
cd /hps/software/users/birney/ian/repos/pilot_paper

singularity build --remote \
    $CONT \
    code/snakemake/20210701/workflow/envs/R_4.1.0/R_4.1.0.def 

ssh proxy-codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
singularity shell --bind /hps/software/users/birney/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/software/users/birney/ian/tmp:/tmp \
                  --bind /hps/software/users/birney/ian/run:/run \
                  $CONT

# Then run rserver, setting path of config file containing library path
rserver --rsession-config-file /hps/software/users/birney/ian/repos/pilot_paper/code/snakemake/20210701/workflow/envs/rstudio_server/rsession.conf

# Then in another terminal
ssh -L 8787:hl-codon-44-04:8787 proxy-codon