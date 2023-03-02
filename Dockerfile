FROM continuumio/miniconda3
LABEL maintainer="lachlan.grose@monash.edu"
#This docker image has been adapted from the lavavu dockerfile
# install things

RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
    gcc \
    g++ \
    libc-dev \
    gfortran \
    openmpi-bin \
    libopenmpi-dev \ 
    make 
# RUN conda install -c conda-forge python=3.9 -y
RUN conda install -c conda-forge -c loop3d\
    pip \
    map2model\
    hjson\
    owslib\
    beartype\
    gdal=3.5.2\
    rasterio=1.2.10 \
    meshio\
    scikit-learn \
    cython \
    numpy \
    pandas \
    scipy \
    pymc3 \
    jupyter \
    pyamg \
    # arviz==0.11.0 \
    pygraphviz \
    geopandas \
    shapely \
    ipywidgets \
    ipyleaflet \
    folium \
    jupyterlab \
    nodejs \
    rasterio\
    geopandas\
    -y

RUN pip install ipyfilechooser
RUN jupyter nbextension enable --py --sys-prefix ipyleaflet
RUN pip install lavavu-osmesa mplstereonet

ENV LD_LIBRARY_PATH=/opt/conda/lib/python3.10/site-packages/lavavu/


ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

USER root
RUN chown -R ${NB_UID} ${HOME}

RUN pip install snakeviz

# Add Tini
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
ENTRYPOINT ["/tini", "--"]

USER ${NB_USER}

RUN mkdir notebooks
RUN git clone https://github.com/Loop3D/map2loop-2.git map2loop
RUN git clone https://github.com/Loop3D/LoopProjectFile.git 
RUN git clone https://github.com/TOMOFAST/Tomofast-x.git
RUN pip install LoopStructural
RUN pip install -e map2loop
RUN pip install -e LoopProjectFile
# WORKDIR Tomofast-x
# RUN make
WORKDIR ${HOME}/notebooks

# RUN pip install -e LoopStructural
CMD ["jupyter", "lab", "--ip='0.0.0.0'", "--NotebookApp.token=''", "--no-browser" ]

EXPOSE 8050
EXPOSE 8080:8090