# Base image
FROM ubuntu:latest
WORKDIR /code
RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y build-essential python3-dev
RUN apt-get install -y python3 python3-pip
#RUN python3 install pip --upgrade
ADD requirements.txt /code
RUN pip3 install -r requirements.txt
ADD setup.py /code
RUN mkdir /code/FME
ADD FME /code/FME
RUN python3 setup.py install 
RUN mkdir /notebooks
RUN jupyter notebook --generate-config
RUN apt-get install -y git pybind11-dev
RUN git clone --recursive http://git.tiker.net/trees/meshpy.git meshpy
WORKDIR /code/meshpy
RUN ls
RUN ./configure
RUN python3 setup.py install
RUN echo "c.NotebookApp.password='loop3d'">>/root/.jupyter/jupyter_notebook_config.p
CMD jupyter notebook --no-browser --ip 0.0.0.0 --port 8888 /notebooks --allow-root c.NotebookApp.token='loop3d'