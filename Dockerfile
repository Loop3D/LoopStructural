# Base image
FROM ubuntu:latest
WORKDIR /code
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential python3-dev
RUN apt-get install -y python3 python3-pip
#RUN python3 install pip --upgrade
# install required generic python packages
RUN pip3 install cython numpy jupyter matplotlib scipy
ADD setup.py /code
RUN mkdir /code/FME
ADD FME /code/FME
RUN python3 setup.py install 
RUN mkdir /notebooks
RUN jupyter notebook --generate-config
RUN apt-get install -y git pybind11-dev 
RUN git clone --recursive http://git.tiker.net/trees/meshpy.git meshpy
WORKDIR /code/meshpy
# install FME dependencies pybind11/meshypy
RUN pip3 install pybind11
RUN apt-get install -y python3-pybind11 libx11-dev libpng-dev
RUN python3 setup.py install
# install useful libraries
RUN apt-get install -y libgl-dev libsuitesparse-dev libsm6
RUN pip3 install geopandas vista scikit-learn meshio pyevtk vtk lavavu scikit-sparse jupyter_contrib_nbextensions
RUN jupyter contrib nbextension install --user
RUN echo "c.NotebookApp.token='loop3d'">>/root/.jupyter/jupyter_notebook_config.py
CMD jupyter notebook --no-browser --ip 0.0.0.0 --port 8888 /notebooks --allow-root c.NotebookApp.token='loop3d'