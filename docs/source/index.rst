.. LoopStructural documentation master file, created by
   sphinx-quickstart on Thu May 28 15:51:42 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

LoopStructural
===============
.. image:: ./images/image823.png


 
Overview
========

LoopStructural is an open source 3D modelling library providing access to multiple interpolation schemes  with a
high level and easy to use API for creating geological models. The library has been written for the Loop platform by
Lachlan Grose at Monash University.

Loop is an open source 3D probabilistic geological and geophysical modelling platform,
initiated by Geoscience Australia and the OneGeology consortium. The project is funded by Australian territory,
State and Federal Geological Surveys, the Australian Research Council and the MinEx Collaborative Research Centre.

Examples
==========

LoopStructural can be used to build a 3D geological model using geological relationships between geological objects
e.g. faults, folds, unconformities and stratigraphic contacts. The library also provides a high level API to access
the fast interpolation schemes that are used by LoopStructural.

Using GeologicalModel
----------------------

The following example shows how to use the geological model interface to create a geological model from a dataset and 
evaluate the scalar field and gradient of the interpolator at some random locations.

.. pyvista-plot::
   :context:
   :include-source: true
   :force_static:

    from LoopStructural import GeologicalModel
    from LoopStructural.datatypes import BoundingBox
    from LoopStructural.visualisation import Loop3DView
    from LoopStructural.datasets import load_claudius
    
    import numpy as np
    data, bb = load_claudius()
    
    #bb constaints origin and maximum of axis aligned bounding box
    #data is a pandas dataframe with X,Y,Z,val,nx,ny,nz, feature_name
    
    model = GeologicalModel(bb[0,:],bb[1,:])
    model.data = data
    # nelements specifies the number of discrete interpolation elements
    # 'stratí' is the feature name in the data dataframe
    model.create_and_add_foliation('strati',nelements=1e5)
    model.update()
    # get the value of the interpolator at some random locations
    locations = np.array(
        [
            np.random.uniform(bb[0, 0], bb[1, 0],5),
            np.random.uniform(bb[0, 1], bb[1, 1],5),
            np.random.uniform(bb[0, 2], bb[1, 2],5),
        ]
    ).T
    val = model.evaluate_feature_value('strati', locations)
    # get the gradient of the interpolator
    gradient = model.evaluate_feature_gradient('strati',locations)
    
    #Plot the scalar field of the model
    model['strati'].scalar_field().plot()


Using InterpolatorBuilder
-------------------------

To access the interpolators directly the `InterpolatorBuilder` can be used
to help assemble an interpolator from a combination of value, gradient and inequality 
constraints.

.. pyvista-plot::
   :context:
   :include-source: true
   :force_static:

    from LoopStructural import BoundingBox, InterpolatorBuilder
    from LoopStructural.utils import EuclideanTransformation
    from LoopStructural.datasets import load_claudius

    # load in a dataframe with x,y,z,val,nx,ny,nz to constrain the value and
    # gradient normal of the interpolator
    data, bb = load_claudius()
    # create a transformer to move the data to the centre of mass of the dataset
    # this avoid any numerical issues caused by large coordinates. This dataset
    # is already axis aligned so we don't need to rotate the data but we can rotate it
    # so that the main anisotropy of the dataset is aligned with the x axis.
    transformer = EuclideanTransformation(dimensions=3, fit_rotation=False).fit(
        data[["X", "Y", "Z"]]
    )
    data[["X", "Y", "Z"]] = transformer.transform(data[["X", "Y", "Z"]])
    # Find the bounding box of the data by finding the extent
    bounding_box = BoundingBox().fit(data[["X", "Y", "Z"]])
    # assemble an interpolator using the bounding box and the data
    builder = (
        InterpolatorBuilder(interpolatortype="FDI", bounding_box=bounding_box)
        .create_interpolator()
        .set_value_constraints(data.loc[data["val"].notna(), ["X", "Y", "Z", "val"]].values)
        .set_normal_constraints(
            data.loc[data["nx"].notna(), ["X", "Y", "Z", "nx", "ny", "nz"]].values
        )
    )
    interpolator = builder.build_interpolator()
    # Set the number of elements in the bounding box to 10000 and create a structured grid
    bounding_box.nelements = 10000
    mesh = bounding_box.structured_grid()
    # add the interpolated values to the mesh at the nodes
    # mesh.properties["v"] = interpolator.evaluate_value(bounding_box.regular_grid(order='F'))
    mesh.properties["v"] = interpolator.evaluate_value(mesh.nodes)

    # or the cell centres
    # mesh.cell_properties["v"] = interpolator.evaluate_value(mesh.cell_centres)


    # visualise the scalar value

    mesh.plot()

    # We can also add gradient properties and visualise these
    mesh.properties["grad"] = interpolator.evaluate_gradient(mesh.nodes)
    mesh.vtk().glyph(orient="grad").plot()

.. toctree::
   :hidden:

   getting_started/index
   _auto_examples/index
   user_guide/index

   
   .. toctree::
   :caption: LoopStructural API
   :hidden:

   API


