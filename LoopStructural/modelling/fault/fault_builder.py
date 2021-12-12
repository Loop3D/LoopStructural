from ..features.structural_frame_builder import StructuralFrameBuilder

import numpy as np
from LoopStructural.utils import getLogger

logger = getLogger(__name__)


class FaultBuilder(StructuralFrameBuilder):
    def __init__(
        self,
        interpolator=None,
        interpolators=None,
        model=None,
        fault_bounding_box_buffer=0.2,
        **kwargs,
    ):
        """A specialised structural frame builder for building a fault

        Parameters
        ----------
        interpolator : GeologicalInterpolator, optional
            the interpolator to use for building the fault frame, by default None
        interpolators : [GeologicalInterpolator, GeologicalInterpolator, GeologicalInterpolator], optional
            a list of interpolators to use for building the fault frame, by default None
        model : GeologicalModel
            reference to the model containing the fault
        fault_bounding_box_buffer: float, default 0.2
            the maximum area around the model domain that a fault is modelled. For high displacement faults this
            may need to be large, smaller values will be result in fewer degrees of freedom = quicker interpolation
        """

        StructuralFrameBuilder.__init__(self, interpolator, interpolators, **kwargs)
        self.model = model
        self.origin = np.array([np.nan, np.nan, np.nan])
        self.maximum = np.array(
            [np.nan, np.nan, np.nan]
        )  # self.model.bounding_box[1, :]
        # define a maximum area to mesh adding buffer to model
        # buffer = .2
        self.minimum_origin = self.model.bounding_box[
            0, :
        ] - fault_bounding_box_buffer * (
            self.model.bounding_box[1, :] - self.model.bounding_box[0, :]
        )
        self.maximum_maximum = self.model.bounding_box[
            1, :
        ] + fault_bounding_box_buffer * (
            self.model.bounding_box[1, :] - self.model.bounding_box[0, :]
        )

    def update_geometry(self, points):
        self.origin = np.nanmin(np.array([np.min(points, axis=0), self.origin]), axis=0)
        self.maximum = np.nanmax(
            np.array([np.max(points, axis=0), self.maximum]), axis=0
        )
        self.origin[self.origin < self.minimum_origin] = self.minimum_origin[
            self.origin < self.minimum_origin
        ]
        self.maximum[self.maximum > self.maximum_maximum] = self.maximum_maximum[
            self.maximum > self.maximum_maximum
        ]

    def create_data_from_geometry(
        self,
        data,
        fault_center,
        normal_vector,
        slip_vector,
        minor_axis=None,
        major_axis=None,
        intermediate_axis=None,
        w=1.0,
        points=True,
    ):
        """Generate the required data for building a fault frame for a fault with the
        specified parameters

        Parameters
        ----------
        data : DataFrame,
            model data
        fault_center : np.array(3)
            x,y,z coordinates of the fault center
        normal_vector : np.array(3)
            x,y,z components of normal vector to fault, single observation usually
            average direction
        slip_vector : np.array(3)
            x,y,z components of slip vector for the fault, single observation usually
            average direction
        minor_axis : double
            distance away from fault for the fault volume
        major_axis : double
            fault extent
        intermediate_axis : double
            fault volume radius in the slip direction
        """
        if major_axis is None:
            fault_trace = data.loc[
                np.logical_and(data["coord"] == 0, data["val"] == 0), ["X", "Y"]
            ].to_numpy()
            distance = np.linalg.norm(
                fault_trace[:, None, :] - fault_trace[None, :, :], axis=2
            )
            if len(distance) == 0 or np.sum(distance) == 0:
                logger.error("There is no fault trace for {}".format(self.name))
                # add any data anyway - usually just orientation data
                self.add_data_from_data_frame(data)
                self.origin = self.model.bounding_box[0, :]
                self.maximum = self.model.bounding_box[1, :]
                return
            major_axis = np.max(distance)
            logger.warning(f"Fault major axis using map length: {major_axis}")

        if minor_axis is None:
            minor_axis = major_axis / 2.0
        if intermediate_axis is None:
            intermediate_axis = major_axis
        normal_vector /= np.linalg.norm(normal_vector)
        slip_vector /= np.linalg.norm(slip_vector)
        # check if slip vector is inside fault plane, if not project onto fault plane
        # if not np.isclose(normal_vector @ slip_vector, 0):
        #     logger.info("{} : projecting slip vector onto fault plane".format(self.name))
        #     slip_vector = np.cross(normal_vector, np.cross(slip_vector ,normal_vector))
        strike_vector = np.cross(normal_vector, slip_vector)
        fault_edges = np.zeros((2, 3))
        fault_tips = np.zeros((2, 3))
        fault_depth = np.zeros((2, 3))
        data.reset_index(inplace=True)
        if fault_center is not None:
            if minor_axis is not None:
                fault_edges[0, :] = fault_center[:3] + normal_vector * minor_axis
                fault_edges[1, :] = fault_center[:3] - normal_vector * minor_axis
                self.update_geometry(fault_edges)

                # choose whether to add points -1,1 to constrain fault frame or a scaled
                # vector
                if points == True:
                    data.loc[
                        len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                    ] = [
                        fault_edges[0, 0],
                        fault_edges[0, 1],
                        fault_edges[0, 2],
                        self.name,
                        1,
                        0,
                        w,
                    ]
                    data.loc[
                        len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                    ] = [
                        fault_edges[1, 0],
                        fault_edges[1, 1],
                        fault_edges[1, 2],
                        self.name,
                        -1,
                        0,
                        w,
                    ]
                if points == False:
                    mask = np.logical_and(data["coord"] == 0, ~np.isnan(data["gx"]))
                    data.loc[mask, ["gx", "gy", "gz"]] /= np.linalg.norm(
                        data.loc[mask, ["gx", "gy", "gz"]], axis=1
                    )[:, None]
                    # scale vector so that the distance between -1 and 1 is the minor axis length
                    data.loc[mask, ["gx", "gy", "gz"]] /= minor_axis * 0.5
            if major_axis is not None:
                fault_tips[0, :] = fault_center[:3] + strike_vector * 0.5 * major_axis
                fault_tips[1, :] = fault_center[:3] - strike_vector * 0.5 * major_axis
                self.update_geometry(fault_tips)
                # we want the tips of the fault to be -1 and 1
                data.loc[
                    len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                ] = [
                    fault_center[0],
                    fault_center[1],
                    fault_center[2],
                    self.name,
                    0,
                    2,
                    w,
                ]
                data.loc[
                    len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                ] = [
                    fault_tips[1, 0],
                    fault_tips[1, 1],
                    fault_tips[1, 2],
                    self.name,
                    -0.5,
                    2,
                    w,
                ]
                data.loc[
                    len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                ] = [
                    fault_tips[0, 0],
                    fault_tips[0, 1],
                    fault_tips[0, 2],
                    self.name,
                    0.5,
                    2,
                    w,
                ]
                strike_vector /= major_axis
            if intermediate_axis is not None:

                fault_depth[0, :] = fault_center[:3] + slip_vector * intermediate_axis
                fault_depth[1, :] = fault_center[:3] - slip_vector * intermediate_axis
                data.loc[
                    len(data), ["X", "Y", "Z", "feature_name", "val", "coord", "w"]
                ] = [
                    fault_center[0],
                    fault_center[1],
                    fault_center[2],
                    self.name,
                    0,
                    1,
                    w,
                ]

                self.update_geometry(fault_depth)
                # TODO need to add data here
                slip_vector /= intermediate_axis
                data.loc[
                    len(data),
                    [
                        "X",
                        "Y",
                        "Z",
                        "feature_name",
                        "nx",
                        "ny",
                        "nz",
                        "val",
                        "coord",
                        "w",
                    ],
                ] = [
                    fault_center[0],
                    fault_center[1],
                    fault_center[2],
                    self.name,
                    slip_vector[0],
                    slip_vector[1],
                    slip_vector[2],
                    0,
                    1,
                    w,
                ]
        self.add_data_from_data_frame(data)
        self.update_geometry(data[["X", "Y", "Z"]].to_numpy())

    def set_mesh_geometry(self, buffer, rotation):
        """set the mesh geometry

        Parameters
        ----------
        buffer : double
            percentage of length to add to edges
        """
        length = np.max(self.maximum - self.origin)
        # for builder in self.builders:
        # all three coordinates share the same support
        self.builders[0].set_interpolation_geometry(
            self.origin - length * buffer, self.maximum + length * buffer, rotation
        )

    def add_splay(self, splay, splayregion=None):
        if splayregion is None:

            def splayregion(xyz):
                pts = (
                    self.builders[0].data[["X", "Y", "Z", "val"]].to_numpy()
                )  # get_value_constraints()
                pts = pts[pts[:, 3] == 0, :]
                # check whether the fault is on the hanging wall or footwall of splay fault

                ext_field = splay[2].evaluate_value(pts[:, :3])
                surf_field = splay[0].evaluate_value(pts[:, :3])
                intersection_value = ext_field[np.nanargmin(np.abs(surf_field))]
                mask = np.zeros(xyz.shape[0], dtype="bool")
                val = splay[2].evaluate_value(xyz)

                if np.nanmedian(ext_field) > intersection_value:
                    mask[~np.isnan(val)] = val[~np.isnan(val)] < intersection_value
                    return mask
                elif np.nanmedian(ext_field) < intersection_value:
                    mask[~np.isnan(val)] = val[~np.isnan(val)] > intersection_value
                    return mask
                else:
                    logger.warning(
                        f"Not adding splay, cannot identify splay overlap region for {self.name} and {splay.name}"
                    )
                    return mask

        self.builders[0].add_equality_constraints(splay, splayregion)
        return splayregion

    def update(self):
        for i in range(3):
            self.builders[i].update()

    def up_to_date(self, callback=None):
        for i in range(3):
            self.builders[i].up_to_date(callback=callback)
