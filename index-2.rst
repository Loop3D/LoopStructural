from LoopStructural import BoundingBox, InterpolatorBuilder
from LoopStructural.utils import EuclideanTransformation
from LoopStructural.datasets import load_claudius

# load in a dataframe with x,y,z,val,nx,ny,nz to constrain the value and
# gradient normal of the interpolator
data, bb = load_claudius()

# Find the bounding box of the data by finding the extent
bounding_box = BoundingBox().fit(data[["X", "Y", "Z"]])
# assemble an interpolator using the bounding box and the data
interpolator = (
    InterpolatorBuilder(
        interpolatortype="FDI", bounding_box=bounding_box.with_buffer(0.1)
    )
    .add_value_constraints(data.loc[data["val"].notna(), ["X", "Y", "Z", "val"]].values)
    .add_normal_constraints(
        data.loc[data["nx"].notna(), ["X", "Y", "Z", "nx", "ny", "nz"]].values
    )
    .build()
)
# Set the number of elements in the bounding box to 10000 and create a structured grid
bounding_box.nelements = 10000
mesh = bounding_box.structured_grid()
# add the interpolated values to the mesh at the nodes
mesh.properties["v"] = interpolator.evaluate_value(mesh.nodes)

# or the cell centres
# mesh.cell_properties["v"] = interpolator.evaluate_value(mesh.cell_centres)


# visualise the scalar value

mesh.plot()

# We can also add gradient properties and visualise these
mesh.properties["grad"] = interpolator.evaluate_gradient(mesh.nodes)
mesh.vtk().glyph(orient="grad").plot()