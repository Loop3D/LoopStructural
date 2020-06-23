Windows using linux subsystem
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
.. container:: toggle

    .. container:: header

        **Installing WSL**

    To setup the windows subsystem for linux you must have administrator rights on your computer.
    Open PowerShell as Administrator (right click and choose run as administrator) and run the following command:

    .. image:: images/powershell_enable_wls.png

    .. code-block:: PowerShell

        Enable-WindowsOptionalFeature -Online -FeatureName Microsoft-Windows-Subsystem-Linux

    Once this command has been executed you may need to restart your computer.

    When the WSL has been enabled you can use the microsoft store to install a linux operating system.
    You can use any operating system you want however, this guide assumes you are using ubuntu 18.04 LTS.

    .. image:: images/ubuntu_microsoft_store.png

    Once you have installed the operating system you will be asked to enter a unix username and password, these do not have
    to be the same as windows - but it is important that you remember them.

    You should then have access to the linux terminal. Which will look something like below:

    .. image:: images/wls_terminal.png

.. container:: toggle

    .. container:: header

        **Installing github for desktop**

    Install a windows github client e.g.:

     * `tortoiseGit` <https://desktop.github.com/>
     * `github dekstop` <https://desktop.github.com/>


    Now clone the LoopStructural repository from the website by clicked clone and open in desktop.

    .. image:: images/githubwebsite_clone.png

    This will then ask you to input a directory for cloning the respository:

    .. image:: images/githubwindows_clone.png

    .. code-block::

        C:/Users/{username}/Documents/Repositories/LoopStructural

    MobaXterm is a terminal client with an X11 server it is best to use this to connect to WLS.
    Download `MobaXterm`<https://mobaxterm.mobatek.net/download.html>

    .. image:: images/moba_xterm.png

.. container:: toggle

    .. container:: header

        **Installing LoopStructural**


    Before you install LoopStructural the ubuntu package manager should be updated and the installed packages should be upgraded.

    Using the linux terminal type in the following commands.
    You can paste into the WLS terminal using


    .. code-block::

        sudo apt-get update && sudo apt-get upgrade

    The dependencies can then be installed:

    * python3
    * python3-dev
    * python3-venv
    * pybind11-dev
    * mesa-common-dev
    * mesa-utils
    * libgl1-mesa-dev
    * g++
    * gcc
    * make

    .. code-block::

        sudo apt-get update  && sudo apt-get install python3 python3-venv python3-dev make pybind11-dev mesa-common-dev mesa-utils libgl1-mesa-dev gcc g++

    It is then recommended to create a new python virtual environment for LoopStructural.

    `Python Virtual Environments: A primer` <https://realpython.com/python-virtual-environments-a-primer/>

    You can create the virtual environment in any location on your computer.

    You can change to the directory where LoopStructural is located by using the following command.


    .. code-block::

        cd /mnt/c/Users/{username}/Documents/Repositories/LoopStructural


    Remember to change the path to the directory where LoopStructural is located on your computer.

    You can then create a virtual environment using the following command.
    This creates a virtual environment called venv inside the LoopStructural repository.
    This folder is automatically ignored by git.


    .. code-block::

        python3 -m venv venv

    You can then create an environment variable for this location by editing your .bashrc file.
    Using VIM or your favourite text editor

    .. code-block::

        nano ~/.bashrc

    Add a line to end end of the file:

    .. code-block::

        export LOOP_ENV=/mnt/c/Users/{username}/Documents/Repository/LoopStructural/venv
        alias LoopStructural='. $LOOP_ENV/bin/activate'

    The second line creates a command line command for switching to the LoopStructural virtual environment.

    .. image:: images/edit_bashrc.png


    For convenience you can symbolic link folders to the home directory for linux.
    This means that the LoopStructural folder will appear in the home directory of your linux user.

    .. code-block::

        ln -s /mnt/c/Users/{username}/Documents/Repository/LoopStructural LoopStructural

    Now change directory to the home folder for linux using the terminal

    .. code-block::

        cd ~

    Now change directory into LoopStructural

    .. code-block::

        cd LoopStructural



    You can now install LoopStructural using the makefile.

    .. code-block::

        make all

    This should run the following commands:

    .. code-block::

        . ${LOOP_ENV}/bin/activate  &&
        pip3 install -r requirements.txt &&
        python3 setup.py install build_ext --inplace;


    A jupyter notebook server can be run from within the LoopStructural folder by running

    .. code-block::

        make notebook

    .. image:: images/run_jupyter.png

    You can then navigate to the jupyter notebook server using your browser.


    .. code-block::

        localhost:8888

    .. image:: images/jupyter_browser.png

    You can now start using LoopStructural.
    Try working through one of the examples/tutorials found in the notebooks directory.

.. container:: toggle

    .. container:: header

        **Upgrading LoopStructural**

    If you have already installed LoopStructural and want to upgrade to the most recent version.

    First pull the most recent version from github.

    Using the WSL change to the LoopStructural directory and run the makefile

    .. code-block::

        cd LoopStructural
        pip uninstall LoopStructural
        make build

    This will not install the requirements.txt and only call the setup.py file for LoopStructural.

.. container:: toggle

    .. container:: header

        **Running LoopStructural**

    To use the model viewing capabilities of LavaVu you need to use MobaXterm or another ssh/terminal client with x forwarding capabilities.
    To run the included examples in LoopStructural you can simply run

    .. code-block::

        make notebook

    and then using your web browser navigate to localhost:8888 or whichever port the jupyter notebook server is on.

    If you want to run a jupyter notebook server from another directory you must first activate the LoopStructural python environment.

    .. code-block::

        LoopStructural

    You can then start a jupyter notebook server

    .. code-block::

        jupyter-notebook --no-browser

    You can then navigate to localhost:8888 or the port specified.


