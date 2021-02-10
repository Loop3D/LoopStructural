FROM continuumio/miniconda3
LABEL maintainer="lachlan.grose@monash.edu"
#This docker image has been adapted from the lavavu dockerfile
# install things
RUN apt-get update -qq 

RUN conda install -c conda-forge pip scikit-learn cython numpy pandas scipy pymc3 jupyter -y
RUN pip install lavavu-osmesa

ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

WORKDIR ${HOME}

RUN mkdir shared_volume
run mkdir LoopStructural
