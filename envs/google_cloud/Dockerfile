FROM tensorflow/tensorflow:1.15.5-gpu-py3-jupyter

MAINTAINER Ian Brettell <brettell@ebi.ac.uk>

ENV export LC_ALL="en_US.UTF-8"
ENV export LC_CTYPE="en_US.UTF-8"
ENV export PATH="/opt/conda/bin:$PATH"
ENV export PATH="/opt/conda/envs/idtrackerai/bin:$PATH"

RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install \
        build-essential \
        python-pip \
        wget \
        bzip2 \
        ca-certificates \
        libglib2.0-0 \
        libxext6 \
        libsm6 \
        libxrender1 \
        git \
        locales
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get clean

# Set locale

RUN locale-gen en_US.UTF-8

# Install idtrackerai 

RUN pip install idtrackerai[gpu]

# Configure for Kubeflow Jupyter notebooks

ENV NB_PREFIX /

CMD ["sh","-c", "jupyter notebook --notebook-dir=/home/jovyan --ip=0.0.0.0 --no-browser --allow-root --port=8888 --NotebookApp.token='' --NotebookApp.password='' --NotebookApp.allow_origin='*' --NotebookApp.base_url=${NB_PREFIX}"]

