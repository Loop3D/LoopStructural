from typing import Union
from ._structural_frame_builder import StructuralFrameBuilder
from LoopStructural.utils import get_vectors
import numpy as np
import pandas as pd
from ....utils import getLogger, BoundingBox

logger = getLogger(__name__)


class FaultBuilder(StructuralFrameBuilder):
    def __init__(
        self,
        interpolatortype: Union[str, list],
        bounding_box: BoundingBox,
        nelements: Union[int, list] = 1000,
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
        from LoopStructural.modelling.features.fault import (
            FaultSegment,
        )  # defer import until needed

        StructuralFrameBuilder.__init__(
            self,
            interpolatortype,
            bounding_box,
            nelements,
            frame=FaultSegment,
            **kwargs,
        )
        self.frame.model = model
        self.model = model
        self.origin = np.array([np.nan, np.nan, np.nan])
        self.maximum = np.array(
            [np.nan, np.nan, np.nan]
        )  # self.model.bounding_box[1, :]
        # define a maximum area to mesh adding buffer to model
        # buffer = .2
        self.minimum_origin = bounding_box.with_buffer(fault_bounding_box_buffer).origin
        self.maximum_maximum = bounding_box.with_buffer(
            fault_bounding_box_buffer
        ).maximum

        self.fault_normal_vector = None
        self.fault_slip_vector = None
        self.fault_strike_vector = None
        self.fault_minor_axis = None
        self.fault_major_axis = None
        self.fault_intermediate_axis = None
        self.fault_centre = None

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
        fault_frame_data: pd.DataFrame,
        fault_center=None,
        fault_normal_vector=None,
        fault_slip_vector=None,
        minor_axis=None,
        major_axis=None,
        intermediate_axis=None,
        w=1.0,
        points=False,
        force_mesh_geometry=False,
        fault_buffer=0.2,
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
        trace_mask = np.logical_and(
            fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0
        )
        logger.info(f"There are {np.sum(trace_mask)} points on the fault trace")
        if np.sum(trace_mask) == 0:
            logger.error(
                "You cannot model a fault without defining the location of the fault"
            )
            raise ValueError(f"There are no points on the fault trace")

        # get all of the gradient data associated with the fault trace
        if fault_normal_vector is None:
            gradient_mask = np.logical_and(
                fault_frame_data["coord"] == 0, ~np.isnan(fault_frame_data["gz"])
            )
            vector_data = fault_frame_data.loc[
                gradient_mask, ["gx", "gy", "gz"]
            ].to_numpy()
            normal_mask = np.logical_and(
                fault_frame_data["coord"] == 0, ~np.isnan(fault_frame_data["nz"])
            )
            vector_data = np.vstack(
                [
                    vector_data,
                    fault_frame_data.loc[normal_mask, ["nx", "ny", "nz"]].to_numpy(),
                ]
            )
            if len(vector_data) == 0:
                logger.error(
                    "You cannot model a fault without defining the orientation of the fault\n\
                    Either add orientation data or define the fault normal vector argument"
                )
                raise ValueError(f"There are no points on the fault trace")
            fault_normal_vector = np.mean(vector_data, axis=0)

        logger.info(f"Fault normal vector: {fault_normal_vector}")

        # estimate the fault slip vector

        if fault_slip_vector is None:
            slip_mask = np.logical_and(
                fault_frame_data["coord"] == 1, ~np.isnan(fault_frame_data["gz"])
            )
            fault_slip_data = fault_frame_data.loc[slip_mask, ["gx", "gy", "gz"]]
            if len(fault_slip_data) == 0:
                logger.warning(
                    "There is no slip vector data for the fault, using vertical slip vector\n\
                          projected onto fault surface estimating from fault normal"
                )
                strike_vector, dip_vector = get_vectors(fault_normal_vector[None, :])
                fault_slip_vector = dip_vector[:, 0]
                logger.info(f"Estimated fault slip vector: {fault_slip_vector}")

            fault_slip_vector = fault_slip_data.mean(axis=0).to_numpy()
        if fault_center is None:
            trace_mask = np.logical_and(
                fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0
            )
            fault_center = (
                fault_frame_data.loc[trace_mask, ["X", "Y", "Z"]]
                .mean(axis=0)
                .to_numpy()
            )
        if np.any(np.isnan(fault_slip_vector)):
            logger.info("Fault slip vector is nan, estimating from fault normal")

        self.fault_normal_vector = fault_normal_vector
        self.fault_slip_vector = fault_slip_vector

        self.fault_centre = fault_center
        if major_axis is None:
            fault_trace = fault_frame_data.loc[
                np.logical_and(
                    fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0
                ),
                ["X", "Y"],
            ].to_numpy()
            distance = np.linalg.norm(
                fault_trace[:, None, :] - fault_trace[None, :, :], axis=2
            )
            if len(distance) == 0 or np.sum(distance) == 0:
                logger.warning("There is no fault trace for {}".format(self.name))
                # this can mean there is only a single data point for
                # the fault, its not critical
                # but probably means the fault isn't well defined.
                # add any data anyway - usually just orientation data
                self.add_data_from_data_frame(fault_frame_data)
                self.origin = self.model.bounding_box.origin
                self.maximum = self.model.bounding_box.maximum
                return
            major_axis = np.max(distance)
            logger.warning(f"Fault major axis using map length: {major_axis}")

        if minor_axis is None:
            logger.info(
                f"Fault minor axis not set, using half major axis: {major_axis/2}"
            )
            minor_axis = major_axis / 2.0
        if intermediate_axis is None:
            intermediate_axis = major_axis
            logger.info(
                f"Fault intermediate axis not set, using major axis: {intermediate_axis}"
            )
        self.fault_minor_axis = minor_axis
        self.fault_major_axis = major_axis
        self.fault_intermediate_axis = intermediate_axis
        fault_normal_vector /= np.linalg.norm(fault_normal_vector)
        fault_slip_vector /= np.linalg.norm(fault_slip_vector)
        # check if slip vector is inside fault plane, if not project onto fault plane
        # if not np.isclose(normal_vector @ slip_vector, 0):

        strike_vector = np.cross(fault_normal_vector, fault_slip_vector)
        self.fault_strike_vector = strike_vector

        fault_edges = np.zeros((2, 3))
        fault_tips = np.zeros((2, 3))
        fault_depth = np.zeros((2, 3))
        fault_frame_data.reset_index(inplace=True)
        if not self.fault_major_axis:
            logger.warning(
                f"Fault major axis is not set and cannot be determined from the fault trace. \
            This will result in a fault that is represented by a 1 unit major axis. \
            If this is not intended add major_axis to fault parameters."
            )
        if not self.fault_intermediate_axis:
            logger.warning(
                f"Fault intermediate axis is not set and cannot be determined from the fault trace. \
            This will result in a fault that is represented by a 1 unit intermediate axis. \
            If this is not intended add intermediate_axis to fault parameters."
            )
        if not self.fault_minor_axis:
            logger.warning(
                f"Fault minor axis is not set and cannot be determined from the fault trace. \
            This will result in a fault that is represented by a 1 unit minor axis. \
            If this is not intended add minor_axis to fault parameters."
            )
        if fault_center is not None:
            if minor_axis is not None:
                fault_edges[0, :] = fault_center[:3] + normal_vector * minor_axis
                fault_edges[1, :] = fault_center[:3] - normal_vector * minor_axis
                self.update_geometry(fault_edges)

                # choose whether to add points -1,1 to constrain fault frame or a scaled
                # vector
                if points:
                    fault_frame_data.loc[
                        len(fault_frame_data),
                        ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
                    ] = [
                        fault_edges[0, 0],
                        fault_edges[0, 1],
                        fault_edges[0, 2],
                        self.name,
                        1,
                        0,
                        w,
                    ]
                    fault_frame_data.loc[
                        len(fault_frame_data),
                        ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
                    ] = [
                        fault_edges[1, 0],
                        fault_edges[1, 1],
                        fault_edges[1, 2],
                        self.name,
                        -1,
                        0,
                        w,
                    ]
                    logger.info("Converting fault norm data to gradient data")
                    mask = np.logical_and(
                        fault_frame_data["coord"] == 0,
                        ~np.isnan(fault_frame_data["nx"]),
                    )
                    fault_frame_data.loc[
                        mask, ["gx", "gy", "gz"]
                    ] = fault_frame_data.loc[mask, ["nx", "ny", "nz"]]

                    fault_frame_data.loc[mask, ["nx", "ny", "nz"]] = np.nan
                    mask = np.logical_and(
                        fault_frame_data["coord"] == 0,
                        ~np.isnan(fault_frame_data["gx"]),
                    )
                    fault_frame_data.loc[mask, ["gx", "gy", "gz"]] /= minor_axis * 0.5
                if points == False:
                    logger.info(
                        "Rescaling fault norm constraint length for fault frame"
                    )
                    mask = np.logical_and(
                        fault_frame_data["coord"] == 0,
                        ~np.isnan(fault_frame_data["gx"]),
                    )
                    fault_frame_data.loc[mask, ["gx", "gy", "gz"]] /= np.linalg.norm(
                        fault_frame_data.loc[mask, ["gx", "gy", "gz"]], axis=1
                    )[:, None]
                    # scale vector so that the distance between -1
                    # and 1 is the minor axis length
                    fault_frame_data.loc[mask, ["gx", "gy", "gz"]] /= minor_axis * 0.5
                    mask = np.logical_and(
                        fault_frame_data["coord"] == 0,
                        ~np.isnan(fault_frame_data["nx"]),
                    )
                    fault_frame_data.loc[mask, ["nx", "ny", "nz"]] /= np.linalg.norm(
                        fault_frame_data.loc[mask, ["nx", "ny", "nz"]], axis=1
                    )[:, None]
                    # scale vector so that the distance between -1
                    # and 1 is the minor axis length
                    fault_frame_data.loc[mask, ["nx", "ny", "nz"]] /= minor_axis * 0.5
                    # self.builders[0].add_orthogonal_feature(self,
                    #  feature, w=1.0, region=None, step=1, B=0):

            if major_axis is not None:
                fault_tips[0, :] = fault_center[:3] + strike_vector * 0.5 * major_axis
                fault_tips[1, :] = fault_center[:3] - strike_vector * 0.5 * major_axis
                self.update_geometry(fault_tips)
                # we want the tips of the fault to be -1 and 1
                fault_frame_data.loc[
                    len(fault_frame_data),
                    ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
                ] = [
                    fault_center[0],
                    fault_center[1],
                    fault_center[2],
                    self.name,
                    0,
                    2,
                    w,
                ]
                fault_frame_data.loc[
                    len(fault_frame_data),
                    ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
                ] = [
                    fault_tips[1, 0],
                    fault_tips[1, 1],
                    fault_tips[1, 2],
                    self.name,
                    -0.5,
                    2,
                    w,
                ]
                fault_frame_data.loc[
                    len(fault_frame_data),
                    ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
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
                fault_frame_data.loc[
                    len(fault_frame_data),
                    ["X", "Y", "Z", "feature_name", "val", "coord", "w"],
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
                fault_frame_data.loc[
                    len(fault_frame_data),
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
        self.add_data_from_data_frame(fault_frame_data)
        if force_mesh_geometry:
            self.set_mesh_geometry(fault_buffer, None)
        self.update_geometry(fault_frame_data[["X", "Y", "Z"]].to_numpy())

    def set_mesh_geometry(self, buffer, rotation):
        """set the mesh geometry

        Parameters
        ----------
        buffer : double
            percentage of length to add to edges
        """
        length = np.max(self.maximum - self.origin)
        print(self.origin, length * buffer, self.maximum)
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

        scalefactor = splay.fault_major_axis / self.fault_major_axis
        self.builders[0].add_equality_constraints(splay, splayregion, scalefactor)
        return splayregion

    def update(self):
        for i in range(3):
            self.builders[i].update()

    def up_to_date(self, callback=None):
        for i in range(3):
            self.builders[i].up_to_date(callback=callback)
