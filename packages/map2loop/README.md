![GitHub Release](https://img.shields.io/github/v/release/loop3d/map2loop)
[![DOI](https://img.shields.io/static/v1?label=DOI&message=10.5194/gmd-14-5063-2021&color=blue)](https://doi.org/10.5194/gmd-14-5063-2021)
![License](https://img.shields.io/github/license/loop3d/map2loop)
![PyPI - Downloads](https://img.shields.io/pypi/dm/map2loop?label=pip%20downloads)
![Conda Downloads](https://img.shields.io/conda/dn/loop3d/map2loop?label=Conda%20downloads)
[![Testing](https://github.com/Loop3D/map2loop/actions/workflows/linting_and_testing.yml/badge.svg)](https://github.com/Loop3D/map2loop/actions/workflows/linting_and_testing.yml)
[![Build and Deploy Documentation](https://github.com/Loop3D/map2loop/actions/workflows/documentation.yml/badge.svg)](https://github.com/Loop3D/map2loop/actions/workflows/documentation.yml)


# Map2Loop 3.2

Generate 3D geological model inputs from geological maps — a high-level implementation and extension of the original map2loop code developed by Prof. Mark Jessell at UWA. To see an example interactive model built with map2loop and LoopStructural, follow this link:

<a href="http://tectonique.net/models/brockman_syncline.html">3D Model from the Hamersley region, Western Australia</a>

## Install

#### Option 1: Install with Anaconda

This is the simplest and recommended installation process, with:

```bash
conda install -c loop3d -c conda-forge map2loop
```

#### Option 2: Install with pip
Installation with pip will require that GDAL is installed on your system prior to map2loop installation. 
This is because GDAL cannot be installed <a href='https://hackernoon.com/hn-images/1*m4cnTYJWM7Rmpsju8dSHmQ.jpeg'>via pip</a> (at least not with one line of code), and the GDAL installation process will vary depending on your OS. 

For more information on installing gdal, see <a href="https://pypi.org/project/GDAL/">GDAL's Pypi</a> page.

Once GDAL is available on your system, map2loop can be installed with:
```bash
pip install map2loop
```

#### Option 3: From source

```bash
git clone https://github.com/Loop3D/map2loop.git

cd map2loop

conda install gdal

conda install -c loop3d -c conda-forge --file dependencies.txt

pip install .
```

#### Option 4: From source & developer mode:
```bash
git clone https://github.com/Loop3D/map2loop.git

cd map2loop

conda install gdal

conda install -c loop3d -c conda-forge --file dependencies.txt

pip install -e .
```

### Documentation

Map2loop's documentation is available <a href="https://loop3d.org/map2loop/">here</a>


## Usage

Our notebooks cover use cases in more detail, but here is an example of processing Loop's South Australia remote geospatial data in just 20 lines of Python.

First, let's import map2loop and define a bounding box. You can use GIS software to find one or use [Loop's Graphical User Interface](https://loop3d.github.io/downloads.html) for the best experience and complete toolset. Remember what projection your coordinates are in!

```python
from map2loop.project import Project
from map2loop.m2l_enums import VerboseLevel

# Note that this region is defined in the EPSG 28354 projection and
# x and y represent easting and northing respectively
bbox_3d = {
    'minx': 250805.1529856466,
    'miny': 6405084.328058686,
    'maxx': 336682.921539395,
    'maxy': 6458336.085975628,
    'base': -3200,
    'top': 1200
}
```

Then, specify: the state, directory for the output, the bounding box and projection from above - and hit go! That's it.

```python
proj = Project(use_australian_state_data = "SA",
               working_projection = 'EPSG:28354',
               bounding_box = bbox_3d,
               loop_project_filename = "output.loop3d"
               )

proj.run_all()
```

This is a minimal example and a small part of Loop.

Our _documentation and other resources outline how to extend map2loop and port to the LoopStructural modelling engine. We are working to incorporate geophysical tools and best provide visualisation and workflow consolidation in the GUI._

_Loop is led by Laurent Ailleres (Monash University) with a team of Work Package leaders from:_

- _Monash University: Roy Thomson, Lachlan Grose and Robin Armit_
- _University of Western Australia: Mark Jessell, Jeremie Giraud, Mark Lindsay and Guillaume Pirot_
- _Geological Survey of Canada: Boyan Brodaric and Eric de Kemp_

---

### Known Issues and FAQs

- Developing with docker on Windows means you won't have GPU passthrough and can’t use a discrete graphics card in the container even if you have one.
- If Jupyter links require a token or password, it may mean port 8888 is already in use. To fix, either make docker map to another port on the host ie -p 8889:8888 or stop any other instances on 8888.

### Links

[https://loop3d.github.io/](https://loop3d.github.io/)

[https://github.com/Loop3D/LoopStructural](https://github.com/Loop3D/LoopStructural)
