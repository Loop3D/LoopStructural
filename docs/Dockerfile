FROM continuumio/miniconda3
LABEL maintainer="lachlan.grose@monash.edu"
#This docker image has been adapted from the lavavu dockerfile
# install things

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    gcc \
    g++ \
    libc-dev \
    make\
    libgl1\
    libtinfo5\
    libtiff6\
    libgl1-mesa-glx\
    xvfb
RUN conda install -c conda-forge\
    -c loop3d\
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
    python=3.10\
    pydata-sphinx-theme\
    pyvista\
    loopstructuralvisualisation\
    -y
RUN pip install git+https://github.com/geopandas/geopandas.git@v0.10.2
RUN pip install geoh5py
RUN pip install sphinxcontrib-bibtex
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

RUN mkdir LoopStructural
