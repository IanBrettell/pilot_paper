# For snakemake

# NOTE: raw data locations:
# FTP: /nfs/ftp/private/birney-res-ftp/upload/medaka/videos/ian_pilot/all

####################
# Codon cluster
####################

ssh codon
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
bsub -M 20000 -Is bash
# If needing to copy videos from FTP (rule copy_videos),
# Need to use the datamover queue so that it can see the FTP drive:
# bsub -M 20000 -q datamover -Is bash
cd /hps/software/users/birney/ian/repos/pilot_paper
conda activate snakemake_6.15.5
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

####################
# Fiji
####################

ssh codon
bsub -M 20000 -q gui -XF -Is bash
/hps/nobackup/birney/users/ian/software/Fiji.app/ImageJ-linux64 


#######################
# Build custom containers
#######################

# R
RCONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/R_4.1.2.sif
singularity build --remote \
    $RCONT \
    workflow/envs/R_4.1.2.def

# Open CV (python)
OPENCVCONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/opencv_4.5.1.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $OPENCVCONT \
    workflow/envs/opencv_4.5.1.def

# idtrackerai
IDCONT=/hps/nobackup/birney/users/ian/containers/MIKK_F0_tracking/idtrackerai.sif
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
singularity build --remote \
    $IDCONT \
    workflow/envs/idtrackerai.def

####################
# Run RStudio Server
####################

ssh proxy-codon
bsub -q datamover -M 50000 -Is bash
module load singularity-3.7.0-gcc-9.3.0-dp5ffrp
RCONT=/hps/nobackup/birney/users/ian/containers/pilot_paper/R_4.1.2.sif
singularity shell --bind /hps/nobackup/birney/users/ian/rstudio_db:/var/lib/rstudio-server \
                  --bind /hps/nobackup/birney/users/ian/tmp:/tmp \
                  --bind /hps/nobackup/birney/users/ian/run:/run \
                  $CONT

rserver \
    --rsession-config-file /hps/software/users/birney/ian/repos/pilot_paper/workflow/envs/rsession.conf \
    --server-user brettell

ssh -L 8787:hl-codon-37-04:8787 proxy-codon

####################
# Copy videos from cluster to local
####################

# To set tracking parameters
rsync -aP brettell@codon:/nfs/research/birney/users/ian/pilot/split ~/Desktop/pilot_videos
