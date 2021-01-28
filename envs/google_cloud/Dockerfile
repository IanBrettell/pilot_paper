FROM ubuntu:18.04

MAINTAINER Ian Brettell <brettell@ebi.ac.uk>

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
      locales
RUN rm -rf /var/lib/apt/lists/*
RUN apt-get clean

# Set locale

RUN locale-gen en_US.UTF-8
