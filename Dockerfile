FROM python:3.7-slim

LABEL maintainer="lachlan.grose@monash.edu"
#This docker image has been adapted from the lavavu dockerfile
# install things
RUN apt-get update -qq && \
    DEBIAN_FRONTEND=noninteractive apt-get install -yq --no-install-recommends \
        bash-completion \
        build-essential \
        xorg-dev \
        ssh \
        curl \
        libfreetype6-dev \
        libpng-dev \
        libtiff-dev \
        libxft-dev \
        xvfb \
        freeglut3 \
        freeglut3-dev \
        libgl1-mesa-dri \
        libgl1-mesa-glx \
        mesa-utils \
        libavcodec-dev \
        libavformat-dev \
        libavutil-dev \
        libswscale-dev \
        rsync \
        xauth \
		git \
		pybind11-dev \
		gcc \
		libgl-dev \
		libsuitesparse-dev \
		libsm6
# Add Tini
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook
COPY requirements.txt ${HOME}
RUN pip install setuptools

RUN pip install -r requirements.txt
RUN pip install jupyter_contrib_nbextensions




ENV NB_USER jovyan
ENV NB_UID 1000
ENV HOME /home/${NB_USER}

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}

# script for xvfb-run.  all docker commands will effectively run under this via the entrypoint
RUN printf "#\041/bin/sh \n rm -f /tmp/.X99-lock && xvfb-run -s '-screen 0 1600x1200x16' \$@" >> /usr/local/bin/xvfbrun.sh && \
chmod +x /usr/local/bin/xvfbrun.sh

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}


RUN jupyter contrib nbextension install --user 
RUN jupyter nbextension enable scratchpad/main 
RUN jupyter nbextension enable comment-uncomment/main  
RUN jupyter nbextension enable collapsible_headings/main  

#Install LoopStructural including cython
# setup environment
ENV PYTHONPATH $PYTHONPATH:${HOME}
USER root
RUN python setup.py install build_ext --inplace
USER ${NB_USER}

# Add a notebook profile.
RUN cd ~ && \
    echo "c.NotebookApp.ip = '0.0.0.0'" >> .jupyter/jupyter_notebook_config.py && \
    echo "c.NotebookApp.token = ''" >> .jupyter/jupyter_notebook_config.py

# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/tini", "--", "/usr/local/bin/xvfbrun.sh"]
USER root
RUN rm -r ${HOME}/*
#COPY ../SGTSG_2019/notebooks ${HOME}
RUN mkdir shared_volume
#COPY 
USER ${NB_USER}
# launch notebook
CMD ["jupyter", "notebook", "--ip='0.0.0.0'", "--NotebookApp.token='' ", "--no-browser" ]
