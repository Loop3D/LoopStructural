from typing import Union

from LoopStructural.utils.maths import rotation
from ._structural_frame_builder import StructuralFrameBuilder
from .. import AnalyticalGeologicalFeature
import numpy as np
import pandas as pd
from ....utils import getLogger
from ....datatypes import BoundingBox

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
            model=model,
            **kwargs,
        )
        self.frame.model = model
        self.origin = np.array([np.nan, np.nan, np.nan])
        self.maximum = np.array([np.nan, np.nan, np.nan])  # self.model.bounding_box[1, :]
        self.raw_data = pd.DataFrame()
        if bounding_box is None:
            raise ValueError("BoundingBox cannot be None")

        # define a maximum area to mesh adding buffer to model
        self.minimum_origin = bounding_box.with_buffer(fault_bounding_box_buffer).origin
        self.maximum_maximum = bounding_box.with_buffer(fault_bounding_box_buffer).maximum

        # Use private attributes and property setters to auto-rebuild when geometrical
        # attributes change. Suspend rebuild during initialization.
        self._suspend_rebuild = True
        self._fault_normal_vector = None
        self._fault_slip_vector = None
        self._fault_strike_vector = None
        self._fault_minor_axis = None
        self._fault_major_axis = None
        self._fault_intermediate_axis = None
        self._fault_centre = None
        self._suspend_rebuild = False

    def update_geometry(self, points):
        """
        Update the geometry of the fault by adjusting the origin and maximum bounds
        based on the provided points.

        Parameters
        ----------
        points : numpy.ndarray
            Array of points used to update the fault geometry.
        """
        if self.origin is None or self.maximum is None:
            raise ValueError("Origin and maximum must be initialized before updating geometry.")

        self.origin = np.nanmin(np.array([np.min(points, axis=0), self.origin]), axis=0)
        self.maximum = np.nanmax(np.array([np.max(points, axis=0), self.maximum]), axis=0)
        # add a small buffer 10% of current length to the origin and maximum
        self.origin = self.origin - 0.1 * (self.maximum - self.origin)
        self.maximum = self.maximum + 0.1 * (self.maximum - self.origin)
        self.origin[self.origin < self.minimum_origin] = self.minimum_origin[
            self.origin < self.minimum_origin
        ]
        self.maximum[self.maximum > self.maximum_maximum] = self.maximum_maximum[
            self.maximum > self.maximum_maximum
        ]

    def _on_geometry_change(
        self,
        force_mesh_geometry: bool = False,
        fault_buffer: float | None = None,
        fault_trace_anisotropy: float = 1.0,
        fault_dip_anisotropy: float = 1.0,
        fault_dip: float | None = None,
        fault_pitch: float | None = None,
    ):
        """Rebuild data from the current stored geometrical attributes.

        This calls create_data_from_geometry with a copy of the last raw_data so that
        changing properties triggers a rebuild.
        """
        if getattr(self, "raw_data", None) is None or self.raw_data.empty:
            return

        # Use a copy of the DataFrame so callers don't observe mutations
        self.create_data_from_geometry(
            self.raw_data.copy(),
            fault_center=self._fault_centre,
            fault_normal_vector=self._fault_normal_vector,
            fault_slip_vector=self._fault_slip_vector,
            minor_axis=self._fault_minor_axis,
            major_axis=self._fault_major_axis,
            intermediate_axis=self._fault_intermediate_axis,
            w=1.0,
            points=False,
            force_mesh_geometry=force_mesh_geometry,
            fault_buffer=fault_buffer if fault_buffer is not None else 0.2,
            fault_trace_anisotropy=fault_trace_anisotropy,
            fault_dip=fault_dip if fault_dip is not None else getattr(self, "fault_dip", 90),
            fault_dip_anisotropy=fault_dip_anisotropy,
            fault_pitch=fault_pitch,
        )

    @property
    def fault_normal_vector(self):
        return self._fault_normal_vector

    @fault_normal_vector.setter
    def fault_normal_vector(self, v):
        if v is None:
            self._fault_normal_vector = None
        else:
            if len(v.shape) != 1 or v.shape[0] != 3:
                raise ValueError("fault_normal_vector must be a 3 element array")
            if v.ndim != 1:
                v = v[0,:]
            if v.shape[0] != 3:
                raise ValueError("fault_normal_vector must be a 3 element array")
            arr = np.array(v, dtype=float)
            norm = np.linalg.norm(arr)
            if norm == 0:
                raise ValueError("fault_normal_vector cannot be the zero vector")
            self._fault_normal_vector = arr
        # trigger rebuild unless suspended
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_slip_vector(self):
        return self._fault_slip_vector

    @fault_slip_vector.setter
    def fault_slip_vector(self, v):
        if v is None:
            self._fault_slip_vector = None
        else:
            arr = np.array(v, dtype=float)
            norm = np.linalg.norm(arr)
            if norm == 0:
                raise ValueError("fault_slip_vector cannot be the zero vector")
            self._fault_slip_vector = arr / norm
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_strike_vector(self):
        return self._fault_strike_vector

    @fault_strike_vector.setter
    def fault_strike_vector(self, v):
        self._fault_strike_vector = None if v is None else np.array(v, dtype=float)
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_minor_axis(self):
        return self._fault_minor_axis

    @fault_minor_axis.setter
    def fault_minor_axis(self, v):
        self._fault_minor_axis = None if v is None else float(v)
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_major_axis(self):
        return self._fault_major_axis

    @fault_major_axis.setter
    def fault_major_axis(self, v):
        self._fault_major_axis = None if v is None else float(v)
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_intermediate_axis(self):
        return self._fault_intermediate_axis

    @fault_intermediate_axis.setter
    def fault_intermediate_axis(self, v):
        self._fault_intermediate_axis = None if v is None else float(v)
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

    @property
    def fault_centre(self):
        return self._fault_centre

    @fault_centre.setter
    def fault_centre(self, v):
        if v is None:
            self._fault_centre = None
        else:
            arr = np.array(v, dtype=float)
            if arr.shape[0] != 3:
                raise ValueError("fault_center must be a 3 element array")
            self._fault_centre = arr
        if not getattr(self, "_suspend_rebuild", False):
            self._on_geometry_change()

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
        fault_trace_anisotropy=1.0,
        fault_dip=90,
        fault_dip_anisotropy=1.0,
        fault_pitch=None,
    ):
        """
        Generate the required data for building a fault frame with the specified parameters.

        Parameters
        ----------
        fault_frame_data : pandas.DataFrame
            DataFrame containing fault frame data.
        fault_center : array-like, optional
            Coordinates of the fault center.
        fault_normal_vector : array-like, optional
            Normal vector of the fault.
        fault_slip_vector : array-like, optional
            Slip vector of the fault.
        minor_axis : float, optional
            Minor axis length of the fault.
        major_axis : float, optional
            Major axis length of the fault.
        intermediate_axis : float, optional
            Intermediate axis length of the fault.
        w : float, default=1.0
            Weighting factor for the fault data.
        points : bool, default=False
            Whether to include points in the fault data.
        force_mesh_geometry : bool, default=False
            Whether to force the use of mesh geometry.
        fault_buffer : float, default=0.2
            Buffer size around the fault.
        fault_trace_anisotropy : float, default=1.0
            Anisotropy factor for the fault trace.
        fault_dip : float, default=90
            Dip angle of the fault in degrees.
        fault_dip_anisotropy : float, default=1.0
            Anisotropy factor for the fault dip.
        fault_pitch : float, optional
            Pitch angle of the fault.
        """
        # Suspend automatic rebuilds while we populate internal attributes to avoid recursion
        prev_suspend = getattr(self, "_suspend_rebuild", False)
        self._suspend_rebuild = True
        try:
            self.raw_data = fault_frame_data.copy()
            fault_trace = fault_frame_data.loc[
                np.logical_and(fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0),
                ["X", "Y"],
            ].to_numpy()
            self.fault_dip = fault_dip
            if fault_normal_vector is None:
                if fault_frame_data.loc[
                np.logical_and(fault_frame_data["coord"] == 0, fault_frame_data["nx"].notna())].shape[0]>0:
                    fault_normal_vector = fault_frame_data.loc[
                        np.logical_and(fault_frame_data["coord"] == 0, fault_frame_data["nx"].notna()),
                        ["nx", "ny", "nz"],
                    ].to_numpy().mean(axis=0)

                else:

                    # Calculate fault strike using eigenvectors
                    pts = fault_trace - fault_trace.mean(axis=0)
                    # Calculate covariance matrix
                    cov_matrix = pts.T @ pts
                    # Get eigenvectors and eigenvalues
                    eigenvalues, eigenvectors = np.linalg.eigh(cov_matrix)
                    # Use eigenvector with largest eigenvalue as strike direction
                    strike_vector = eigenvectors[:, np.argmax(eigenvalues)]
                    strike_vector = np.append(strike_vector, 0)  # Add z component
                    strike_vector /= np.linalg.norm(strike_vector)

                    fault_normal_vector = np.cross(strike_vector, [0, 0, 1])
                    # Rotate the fault normal vector according to the fault dip
                    rotation_matrix = rotation(strike_vector[None, :], np.array([90 - fault_dip]))
                    fault_normal_vector = np.einsum("ijk,ik->ij", rotation_matrix, fault_normal_vector[None, :])[0]

            if not isinstance(fault_normal_vector, np.ndarray):
                fault_normal_vector = np.array(fault_normal_vector)

            if fault_pitch is not None:
                rotation_matrix = rotation(fault_normal_vector[None, :], np.array([fault_pitch]))
                fault_slip_vector = np.einsum("ijk,ik->ij", rotation_matrix, fault_normal_vector[None, :])[0]

            if fault_slip_vector is None:
                if fault_frame_data.loc[
                np.logical_and(fault_frame_data["coord"] == 1, fault_frame_data["nx"].notna())].shape[0]>0:
                    fault_slip_vector = fault_frame_data.loc[
                        np.logical_and(fault_frame_data["coord"] == 1, fault_frame_data["nx"].notna()),
                        ["nx", "ny", "nz"],
                    ].to_numpy().mean(axis=0)

                else:
                    fault_slip_vector = np.cross(fault_normal_vector, [1., 0., 0.])
                    if np.linalg.norm(fault_slip_vector) == 0:
                        fault_slip_vector = np.cross(fault_normal_vector, [0., 1., 0.])
                    fault_slip_vector /= np.linalg.norm(fault_slip_vector)
            if fault_center is None:
                fault_trace = fault_frame_data.loc[
                np.logical_and(fault_frame_data["coord"] == 0, fault_frame_data["val"] == 0),
                ["X", "Y"],
            ].to_numpy()
                fault_center = fault_trace.mean(axis=0)
                fault_center = np.array([fault_center[0], fault_center[1], 0.0])
            if not isinstance(fault_center, np.ndarray):
                fault_center = np.array(fault_center)
            if fault_center.shape[0] != 3:
                raise ValueError("fault_center must be a 3 element array")
            # Assign to properties (these won't trigger _on_geometry_change because suspend flag is True)
            self.fault_normal_vector = fault_normal_vector / np.linalg.norm(fault_normal_vector)
            self.fault_slip_vector = fault_slip_vector / np.linalg.norm(fault_slip_vector)

            self.fault_centre = fault_center
            if major_axis is None:

                distance = np.linalg.norm(fault_trace[:, None, :] - fault_trace[None, :, :], axis=2)
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
                logger.info(f"Fault major axis using map length: {major_axis}")

            if minor_axis is None:
                logger.info(f"Fault minor axis not set, using half major axis: {major_axis/2}")
                minor_axis = major_axis / 2.0
            if intermediate_axis is None:
                intermediate_axis = major_axis
                logger.info(f"Fault intermediate axis not set, using major axis: {intermediate_axis}")
            self.fault_minor_axis = minor_axis
            self.fault_major_axis = major_axis
            self.fault_intermediate_axis = intermediate_axis
            self.fault_normal_vector /= np.linalg.norm(self.fault_normal_vector)
            fault_slip_vector /= np.linalg.norm(fault_slip_vector)
            # check if slip vector is inside fault plane, if not project onto fault plane
            # if not np.isclose(normal_vector @ slip_vector, 0):
            strike_vector = np.cross(self.fault_normal_vector, fault_slip_vector)
            self.fault_strike_vector = strike_vector

            fault_edges = np.zeros((2, 3))
            fault_tips = np.zeros((2, 3))
            fault_depth = np.zeros((2, 3))
            fault_frame_data.reset_index(inplace=True)
            if not self.fault_major_axis:
                logger.warning(
                    "Fault major axis is not set and cannot be determined from the fault trace. \
                This will result in a fault that is represented by a 1 unit major axis. \
                If this is not intended add major_axis to fault parameters."
                )
            if not self.fault_intermediate_axis:
                logger.warning(
                    "Fault intermediate axis is not set and cannot be determined from the fault trace. \
                This will result in a fault that is represented by a 1 unit intermediate axis. \
                If this is not intended add intermediate_axis to fault parameters."
                )
            if not self.fault_minor_axis:
                logger.warning(
                    "Fault minor axis is not set and cannot be determined from the fault trace. \
                This will result in a fault that is represented by a 1 unit minor axis. \
                If this is not intended add minor_axis to fault parameters."
                )
            if fault_center is not None:
                if minor_axis is not None:
                    fault_edges[0, :] = fault_center[:3] + self.fault_normal_vector * minor_axis
                    fault_edges[1, :] = fault_center[:3] - self.fault_normal_vector * minor_axis
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
                        fault_frame_data.loc[mask, ["gx", "gy", "gz"]] = fault_frame_data.loc[
                            mask, ["nx", "ny", "nz"]
                        ]

                        fault_frame_data.loc[mask, ["nx", "ny", "nz"]] = np.nan
                        mask = np.logical_and(
                            fault_frame_data["coord"] == 0,
                            ~np.isnan(fault_frame_data["gx"]),
                        )
                        fault_frame_data.loc[mask, ["gx", "gy", "gz"]] /= minor_axis * 0.5
                    if not points:
                        logger.info("Rescaling fault norm constraint length for fault frame")
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
                        if np.sum(mask) == 0:
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
                                self.fault_normal_vector[0] / minor_axis * 0.5,
                                self.fault_normal_vector[1] / minor_axis * 0.5,
                                self.fault_normal_vector[2] / minor_axis * 0.5,
                                np.nan,
                                0,
                                w,
                            ]

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
                    fault_depth[0, :] = fault_center[:3] + fault_slip_vector * intermediate_axis
                    fault_depth[1, :] = fault_center[:3] - fault_slip_vector * intermediate_axis
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
                    fault_slip_vector /= intermediate_axis
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
                        fault_slip_vector[0],
                        fault_slip_vector[1],
                        fault_slip_vector[2],
                        0,
                        1,
                        w,
                    ]

            self.add_data_from_data_frame(fault_frame_data)
            if fault_trace_anisotropy > 0:
                self.add_fault_trace_anisotropy(w=fault_trace_anisotropy)
            if fault_dip_anisotropy > 0:
                self.add_fault_dip_anisotropy(w=fault_dip_anisotropy)
            if force_mesh_geometry:
                self.origin = self.model.bounding_box.origin
                self.maximum = self.model.bounding_box.maximum
            else:
                self.update_geometry(fault_frame_data[["X", "Y", "Z"]].to_numpy())
            self.set_mesh_geometry(fault_buffer, None)
        finally:
            # restore previous suspend state so property setters behave as before
            self._suspend_rebuild = prev_suspend

    def set_mesh_geometry(self, buffer, rotation):
        """set the mesh geometry

        Parameters
        ----------
        buffer : double
            percentage of length to add to edges
        """
        length = np.nanmax(self.maximum - self.origin)
        # for builder in self.builders:
        # all three coordinates share the same support
        self.builders[0].set_interpolation_geometry(
            self.origin - length * buffer, self.maximum + length * buffer, rotation
        )
        self.builders[1].set_interpolation_geometry(
            self.origin - length * buffer, self.maximum + length * buffer, rotation
        )
        self.builders[2].set_interpolation_geometry(
            self.origin - length * buffer, self.maximum + length * buffer, rotation
        )

    def add_splay(self, splay, splay_region=None):
        if splay_region is None:

            def default_splay_region(xyz):
                pts = (
                    self.builders[0].data["X", "Y", "Z", "val"].to_numpy()
                )  # get_value_constraints()
                pts = pts[pts[:, 3] == 0, :]
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

            splay_region = default_splay_region

        scalefactor = splay.fault_major_axis / self.fault_major_axis
        self.builders[0].add_equality_constraints(splay, splay_region, scalefactor)
        return splay_region

    def add_fault_trace_anisotropy(self, w: float = 1.0):
        """
        Add fault trace anisotropy to the model.

        Parameters
        ----------
        w : float, optional
            Weighting factor for the anisotropy, by default 1.0
        """
        if w > 0:
            if self.fault_normal_vector is None:
                raise ValueError("fault_normal_vector must be initialized before adding anisotropy.")

            plane = np.array([0, 0, 1])
            strike_vector = (
                self.fault_normal_vector - np.dot(self.fault_normal_vector, plane) * plane
            )
            strike_vector /= np.linalg.norm(strike_vector)
            strike_vector = np.array([strike_vector[1], -strike_vector[0], 0])

            anisotropy_feature = AnalyticalGeologicalFeature(
                vector=strike_vector, origin=np.array([0, 0, 0]), name="fault_trace_anisotropy"
            )
            self.builders[0].add_orthogonal_feature(
                anisotropy_feature, w=w, region=None, step=1, B=0
            )

    def add_fault_dip_anisotropy(self, w: float = 1.0):
        """
        Add fault dip anisotropy to the model.

        Parameters
        ----------
        w : float, optional
            Weighting factor for the anisotropy, by default 1.0
        """
        if w > 0:
            if self.fault_normal_vector is None:
                raise ValueError("fault_normal_vector must be initialized before adding anisotropy.")

            plane = np.array([0, 0, 1])
            strike_vector = (
                self.fault_normal_vector - np.dot(self.fault_normal_vector, plane) * plane
            )
            strike_vector /= np.linalg.norm(strike_vector)
            strike_vector = np.array([strike_vector[1], -strike_vector[0], 0])

            dip_vector = np.cross(strike_vector, self.fault_normal_vector)
            dip_vector /= np.linalg.norm(dip_vector)

            anisotropy_feature = AnalyticalGeologicalFeature(
                vector=dip_vector, origin=np.array([0, 0, 0]), name="fault_dip_anisotropy"
            )
            self.builders[0].add_orthogonal_feature(
                anisotropy_feature, w=w, region=None, step=1, B=0
            )

    def update(self):
        for i in range(3):
            self.builders[i].update()

    def up_to_date(self, callback=None):
        for i in range(3):
            self.builders[i].up_to_date(callback=callback)
