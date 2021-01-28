FROM julia:1.5.3

MAINTAINER Ian Brettell <brettell@ebi.ac.uk>

COPY envs/julia/Manifest.toml
COPY envs/julia/Project.toml

RUN apt-get update && apt-get -y upgrade
RUN apt-get -y install \
      build-essential \
      wget \
      bzip2 \
      ca-certificates \
      libglib2.0-0 \
      libxext6 \
      libsm6 \
      libxrender1 \
      git \
      locales \
      curl
