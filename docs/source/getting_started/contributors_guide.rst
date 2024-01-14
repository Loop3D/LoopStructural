Contribute
===============================
Contributions are welcome to LoopStructural. 
Contributions can be bug reports, documentation, examples, questions, feature requests, new features and new code.
Please follow the following guidelines.

Issues
-------
Questions
~~~~~~~~~~
If you have a question about how to use LoopStructural please 
`search <https://github.com/Loop3D/LoopStructural/issues>`_ existing issues and look through the example documentation.
If your question has not been previously addressed then please `submit a new issue <https://github.com/Loop3D/LoopStructural/issues/new/choose>`_.
When submitting a new issue please be as verbose as possible using the provided templates and provide as much information to help us answer your question.

Bug reports
~~~~~~~~~~~
To submit a bug report please use the bug reports template when `submitting an issue <https://github.com/Loop3D/LoopStructural/issues/new/choose>`_.
Please provide us with a `minimal working example <https://en.wikipedia.org/wiki/Minimal_working_example>`_ that allows for the bug
to be easily reproduced this will help us in determining the source of the bug. 

Feature requests
~~~~~~~~~~~~~~~~
If you have an idea of a new feature for LoopStructural please submit this using the `Feature request <https://github.com/Loop3D/LoopStructural/issues/new/choose>`_ template.

Contributing new code
----------------------
Any contributions to the documentation and code for LoopStructural are welcome.

If you would like to contribute code to LoopStructural please open an `issue <https://github.com/Loop3D/LoopStructural/issues/new/choose>`_ either 
a bug report or a feature request and then submit a pull request and link it to this issue.

Getting setup
~~~~~~~~~~~
To get started, fork and clone the repository. It's then best to get started with conda::

    conda activate
    pip install -r requirements.txt
    conda list

N.B. On Linux, the LavaVu package requires the installation of extra tools::

    sudo apt install build-essential libgl1-mesa-dev libx11-dev zlib1g-dev

For changes to ``./docs``::

    pip install -r docs/requirements.txt
    make -C ./docs html

Building the docs the first time takes a little while. After this, you will see a build directory. From here you can preview the generated HTML 
in your browser. 

Commit messages
---------------
A changelog is automatically generated from the commit message, the following keywords are identified in the commit log.
 
 * add - defining a new feature that is added to LoopStructural
 * fix - a bugfix
 * change - any changes to features 
 * remove - removing features or code from LoopStructural
 * merge 
 * doc - changes to documentation

The commit messages should include a short summary in the first line that will be referenced in the change log.
 
Version numbering
-----------------
LoopStructural is versions according to semantic versioning:

Given a version number MAJOR.MINOR.PATCH, increment the:

MAJOR version when you make incompatible API changes,
MINOR version when you add functionality in a backwards compatible manner, and
PATCH version when you make backwards compatible bug fixes.
Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.


Documentation
~~~~~~~~~~~~~~
LoopStructural documentation is produced using sphinx with examples being provided by sphinx gallery. 
The code documentation is written using the numpy docstring format for example a function should be documented using the 
following style:

.. code-block:: default


    def new_function(a, b, c=None,d=None):
        """[summary]

        [extended_summary]

        Parameters
        ----------
        a : [type]
            [description]
        b : [type]
            [description]
        c : [type], optional
            [description], by default None
        d : [type], optional
            [description], by default None

        Returns
        -------
        results : [type]
            [description]
        """

License
~~~~~~~
LoopStructural is licensed unded an MIT license and all contributions MUST conform to this license. 


