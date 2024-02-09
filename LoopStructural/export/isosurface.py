# import numpy as np
# from skimage.measure import marching_cubes


# def extract_isosurface(geologicalfeature, model, isovalue):
#     xyz = model.regular_grid(shuffle=False)
#     vals = model.evaluate_model(xyz, scale=False)

#     verts, faces, normals, values = marching_cubes(
#         vals.reshape(model.nsteps, order='C'), isovalue, spacing=model.step_vector
#     )
#     # verts += np.array([self.bounding_box[0, 0], self.bounding_box[0, 1], self.bounding_box[1, 2]])
#     self.model.rescale(verts)
#     return verts, faces
