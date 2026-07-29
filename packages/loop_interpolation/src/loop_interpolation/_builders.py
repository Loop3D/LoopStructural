# from LoopStructural.utils.exceptions import LoopException
# import numpy as np
# from typing import Optional
# from LoopStructural.interpolators import (
#     P1Interpolator,
#     P2Interpolator,
#     FiniteDifferenceInterpolator,
#     GeologicalInterpolator,
#     DiscreteFoldInterpolator,
#     StructuredGrid,
#     TetMesh,
# )
# from LoopStructural.datatypes import BoundingBox
# from LoopStructural.utils.logging import getLogger

# logger = getLogger(__name__)


# def get_interpolator(
#     bounding_box: BoundingBox,
#     interpolatortype: str,
#     nelements: int,
#     element_volume: Optional[float] = None,
#     buffer: float = 0.2,
#     dimensions: int = 3,
#     support=None,
# ) -> GeologicalInterpolator:
#     # add a buffer to the interpolation domain, this is necessary for
#     # faults but also generally a good
#     # idea to avoid boundary problems
#     # buffer = bb[1, :]
#     origin = bounding_box.with_buffer(buffer).origin
#     maximum = bounding_box.with_buffer(buffer).maximum
#     box_vol = np.prod(maximum - origin)
#     if interpolatortype == "PLI":
#         if support is None:
#             if element_volume is None:
#                 # nelements /= 5
#                 element_volume = box_vol / nelements
#             # calculate the step vector of a regular cube
#             step_vector = np.zeros(3)
#             step_vector[:] = element_volume ** (1.0 / 3.0)
#             # step_vector /= np.array([1,1,2])
#             # number of steps is the length of the box / step vector
#             nsteps = np.ceil((maximum - origin) / step_vector).astype(int)
#             if np.any(np.less(nsteps, 3)):
#                 axis_labels = ["x", "y", "z"]
#                 for i in range(3):
#                     if nsteps[i] < 3:
#                         nsteps[i] = 3
#                         logger.error(
#                             f"Number of steps in direction {axis_labels[i]} is too small, try increasing nelements"
#                         )
#                 logger.error("Cannot create interpolator: number of steps is too small")
#                 raise ValueError("Number of steps too small cannot create interpolator")

#             support = TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)
#         logger.info(
#             "Creating regular tetrahedron mesh with %i elements \n"
#             "for modelling using PLI" % (support.ntetra)
#         )

#         return P1Interpolator(support)
#     if interpolatortype == "P2":
#         if support is not None:
#             logger.info(
#                 "Creating regular tetrahedron mesh with %i elements \n"
#                 "for modelling using P2" % (support.ntetra)
#             )
#             return P2Interpolator(support)
#         else:
#             raise ValueError("Cannot create P2 interpolator without support, try using PLI")

#     if interpolatortype == "FDI":
#         # find the volume of one element
#         if element_volume is None:
#             element_volume = box_vol / nelements
#         # calculate the step vector of a regular cube
#         step_vector = np.zeros(3)
#         step_vector[:] = element_volume ** (1.0 / 3.0)
#         # number of steps is the length of the box / step vector
#         nsteps = np.ceil((maximum - origin) / step_vector).astype(int)
#         if np.any(np.less(nsteps, 3)):
#             logger.error("Cannot create interpolator: number of steps is too small")
#             axis_labels = ["x", "y", "z"]
#             for i in range(3):
#                 if nsteps[i] < 3:
#                     nsteps[i] = 3
#             #         logger.error(
#             #             f"Number of steps in direction {axis_labels[i]} is too small, try increasing nelements"
#             #         )
#             # raise ValueError("Number of steps too small cannot create interpolator")
#         # create a structured grid using the origin and number of steps

#         grid = StructuredGrid(origin=origin, nsteps=nsteps, step_vector=step_vector)
#         logger.info(
#             f"Creating regular grid with {grid.n_elements} elements \n" "for modelling using FDI"
#         )
#         return FiniteDifferenceInterpolator(grid)
#     if interpolatortype == "DFI":
#         if element_volume is None:
#             nelements /= 5
#             element_volume = box_vol / nelements
#         # calculate the step vector of a regular cube
#         step_vector = np.zeros(3)
#         step_vector[:] = element_volume ** (1.0 / 3.0)
#         # number of steps is the length of the box / step vector
#         nsteps = np.ceil((maximum - origin) / step_vector).astype(int)
#         # create a structured grid using the origin and number of steps

#         mesh = TetMesh(origin=origin, nsteps=nsteps, step_vector=step_vector)
#         logger.info(
#             f"Creating regular tetrahedron mesh with {mesh.ntetra} elements \n"
#             "for modelling using DFI"
#         )
#         return DiscreteFoldInterpolator(mesh, None)
#     raise LoopException("No interpolator")
#     # fi interpolatortype == "DFI" and dfi is True:
#     #     if element_volume is None:
#     #         nelements /= 5
#     #         element_volume = box_vol / nelements
#     #     # calculate the step vector of a regular cube
#     #     step_vector = np.zeros(3)
#     #     step_vector[:] = element_volume ** (1.0 / 3.0)
#     #     # number of steps is the length of the box / step vector
#     #     nsteps = np.ceil((bb[1, :] - bb[0, :]) / step_vector).astype(int)
#     #     # create a structured grid using the origin and number of steps
#     #     if "meshbuilder" in kwargs:
#     #         mesh = kwargs["meshbuilder"].build(bb, nelements)
#     #     else:
#     #         mesh = kwargs.get(
#     #             "mesh",
#     #             TetMesh(origin=bb[0, :], nsteps=nsteps, step_vector=step_vector),
#     #         )
#     #     logger.info(
#     #         f"Creating regular tetrahedron mesh with {mesh.ntetra} elements \n"
#     #         "for modelling using DFI"
#     #     )
#     #     return DFI(mesh, kwargs["fold"])
#     # if interpolatortype == "Surfe" or interpolatortype == "surfe":
#     #     # move import of surfe to where we actually try and use it
#     #     if not surfe:
#     #         logger.warning("Cannot import Surfe, try another interpolator")
#     #         raise ImportError("Cannot import surfepy, try pip install surfe")
#     #     method = kwargs.get("method", "single_surface")
#     #     logger.info("Using surfe interpolator")
#     #     return SurfeRBFInterpolator(method)
#     # logger.warning("No interpolator")
#     # raise InterpolatorError("Could not create interpolator")
