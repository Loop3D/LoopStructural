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
    # python"<=3.8"\
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
    # geopandas\
    libstdcxx-ng\
    meshio\
    -y
RUN pip install git+https://github.com/geopandas/geopandas.git@v0.10.2
RUN pip install lavavu-osmesa 
ENV LD_LIBRARY_PATH=/opt/conda/lib/python3.10/site-packages/lavavu_osmesa.libs
RUN conda install -c conda-forge pydata-sphinx-theme 
RUN mkdir LoopStructural
