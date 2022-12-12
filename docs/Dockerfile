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
    cython\
    numpy\ 
    pandas\ 
    imageio\ 
    scipy\ 
    matplotlib\  
    sphinx\
    sphinx-gallery\ 
    sphinxcontrib-bibtex\ 
    sphinx_rtd_theme\ 
    myst-parser\ 
    scikit-learn\ 
    scikit-image\ 
    pyamg\ 
    flake8\
    pytest\
    networkx\
    geopandas\
    meshio\
    -y
RUN pip install lavavu-osmesa 
RUN conda install -c conda-forge pydata-sphinx-theme 
RUN mkdir LoopStructural
