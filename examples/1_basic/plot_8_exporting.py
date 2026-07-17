"""
1h. Exporting models
===============================
Once a model has been built, its surfaces and volumes typically need to be
brought into other software (a GIS, a 3D viewer, another modelling
package). This example shows how to extract and export individual
surfaces from a geological feature.

Supported file formats depend on what is being exported and include
``vtk``, ``ts``/``gocad``, ``obj``, ``json``, ``omf`` and ``geoh5`` -
:code:`save` picks the writer to use from the file extension. The
``geoh5`` and ``omf`` formats are container formats that can store the
whole model (surfaces, block model and data) in a single file; the
geoh5 writer additionally requires the optional ``geoh5py`` package.
"""

import tempfile
import pathlib

from LoopStructural import GeologicalModel
from LoopStructural.datasets import load_claudius

data, bb = load_claudius()

model = GeologicalModel(bb[0, :], bb[1, :])
model.data = data
model.create_and_add_foliation("strati")

# write outputs to a temporary directory so this example doesn't leave
# files behind - replace `output_dir` with a real path to keep the output
output_dir = pathlib.Path(tempfile.mkdtemp())

######################################################################
# Export a single surface
# ~~~~~~~~~~~~~~~~~~~~~~~~
# Isosurfaces can be extracted from a geological feature by calling the
# ``.surfaces()`` method on the feature. The argument is the value, list of
# values, or number of evenly-spaced surfaces to extract. This returns a
# list of :class:`LoopStructural.datatypes.Surface` objects, which expose
# the triangles/vertices/normals directly and can also be written to disk
# with ``.save()``.

surfaces = model['strati'].surfaces(value=0.0)
print(f"{len(surfaces)} surface(s), {len(surfaces[0].vertices)} vertices, "
      f"{len(surfaces[0].triangles)} triangles")

surfaces[0].save(str(output_dir / 'strati_surface.vtk'))

######################################################################
# Exporting multiple horizons
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Passing a list of values (or a count) to ``.surfaces()`` extracts an
# isosurface per value, which can be saved individually - useful for
# exporting each stratigraphic horizon as a separate object.

for i, surface in enumerate(model['strati'].surfaces(value=[0.0, 100.0, 200.0])):
    surface.save(str(output_dir / f'strati_horizon_{i}.vtk'))

print(sorted(p.name for p in output_dir.glob('*')))

######################################################################
# Exporting an entire model
# ~~~~~~~~~~~~~~~~~~~~~~~~~
# ``model.save(filename)`` is intended to walk every stratigraphic and
# fault surface together with the block model and input data, and write
# them all out in one call - one file per object for formats like
# ``vtk``, or everything bundled into a single file for container formats
# like ``geoh5``/``omf``.
#
# .. code:: python
#
#    model.save("model.geoh5")
