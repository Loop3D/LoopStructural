FROM continuumio/miniconda3
LABEL maintainer="lachlan.grose@monash.edu"
#This docker image has been adapted from the lavavu dockerfile
# install things

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    gcc \
    g++ \
    libc-dev \
    make
RUN conda install -c conda-forge\
    python"<=3.9"\
    cython\
    numpy\ 
    pandas\ 
    scipy\ 
    matplotlib\  
    sphinx\
    sphinx-gallery\ 
    myst-parser\ 
    scikit-learn\ 
    scikit-image\ 
    networkx\
    gcc\
    # geopandas\
    # meshio\
    -y
RUN pip install lavavu-osmesa==1.8.32
RUN conda install -c conda-forge pydata-sphinx-theme 
RUN mkdir LoopStructural
