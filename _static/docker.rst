Docker
~~~~~~

.. container:: toggle

    .. container:: header

        **Installing Docker**
    Follow the installation instructions for docker `here <https://docs.docker.com/engine/install/>`_.

.. container:: toggle

    .. container:: header

        **Building and running a docker image**
    Using a github client (e.g. for windows Github Desktop)
    Clone this repository to your local drive and change
    directory to the location where you cloned it using
    :code:`cd ${LOOPSTRUCTURAL_FOLDER}`. The docker
    container can be built by running the following command
    .. code-block::console

        docker build -t=loop .

    LoopStructural can be used by running

    .. code-block::console

        run -i -t -p 8888:8888 loop

    This will start a jupyter notebook server running on :code:`localhost:8888`
    without password or certificate required. Be aware any changes made
    to the notebooks within the notebooks directory will NOT be saved on
    your local versions.

    If you want to use your own data with the docker container you will need
    to link your local directory (this can be anywhere) with the docker container.
    To do this add :code:`-v LOCALDIRPATH:/home/joyvan/shared_volume` to the docker command
    so it becomes :code:`docker run -i -t -p 8888:8888 -v LOCALDIRPATH:/home/joyvan/shared_volume`.
    :code:`LOCALDIRPATH` is the full path to the directory you want to share.