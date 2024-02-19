from skimage.measure import marching_cubes
import pandas as pd
import numpy as np
from ..modelling.features import FeatureType

from ..utils import getLogger

logger = getLogger(__name__)


def calculate_fault_intersections(model, threshold=0.001):
    """Calculate the intersections between faults in the model
    threshold defines the minimum displacement value that is
    considered in a fault volume

    Parameters
    ----------
    model : GeologicalModel
        the model to process
    threshold : float, optional
        displacement threshold for defining fault volume, by default 0.001

    Returns
    -------
    intersections
        dataframe matrix of fault intersections
    """
    fault_names = []
    for f in model.features:
        if f.type == FeatureType.FAULT:
            fault_names.append(f.name)
    fault_matrix = pd.DataFrame(columns=fault_names)
    for f in fault_names:
        fault_matrix.loc[f, :] = 0

    for name in fault_names:
        xyz = model.regular_grid(shuffle=False)
        vals = model[name].evaluate_value(xyz)
        model.step_vector = (model.bounding_box[1, :] - model.bounding_box[0, :]) / model.nsteps
        try:

            verts, faces, normals, values = marching_cubes(
                vals.reshape(model.nsteps, order="C"), 0, spacing=model.step_vector
            )
            verts = model.rescale(verts)
        except (ValueError, RuntimeError) as e:
            print(e)
            logger.warning("Cannot isosurface {} at {}, skipping".format(name, 0))
            continue
        for name2 in fault_names:
            if name2 == name:
                continue
            val = model[name2].evaluate_value(model.scale(verts, inplace=False))

            mask = model[name2].inside_volume(model.scale(verts, inplace=False), threshold)
            val[~mask] = np.nan
            if np.all(np.isnan(val)):
                continue
            if np.all(val[~np.isnan(val)] > 0):
                continue  # print('all > 0')
            elif np.all(val[~np.isnan(val)] < 0):
                continue
            else:
                fault_matrix.loc[name, name2] = 1
    return fault_matrix
