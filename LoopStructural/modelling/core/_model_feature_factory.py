"""Feature-construction orchestration for GeologicalModel (see API.md).

Extracted from GeologicalModel to separate the mechanics of building each
feature type (foliation, fold frame, folded foliation, folded fold frame,
intrusion, domain fault, fault) from feature-container state and the public
``create_and_add_*`` API. GeologicalModel's @public_api ``create_and_add_*``
methods stay defined directly on the class -- their __qualname__ is part of
the CI-checked stable API surface -- and dispatch through
``GeologicalModel.create_and_add_feature`` ->
``FeatureBuilderRegistry.create`` -> the staticmethods here (wired up at the
bottom of ``geological_model.py``).
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Callable

import numpy as np
import pandas as pd

from LoopStructural import LoopStructuralConfig

from ...modelling.features import FeatureType, GeologicalFeature, UnconformityFeature
from ...modelling.features.builders import (
    FaultBuilder,
    FoldedFeatureBuilder,
    GeologicalFeatureBuilder,
    StructuralFrameBuilder,
)
from ...modelling.features.fault import FaultSegment
from ...modelling.features.fold import FoldEvent, FoldFrame
from ...modelling.intrusions import IntrusionBuilder, IntrusionFeature, IntrusionFrameBuilder
from ...utils import getLogger

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from .geological_model import GeologicalModel

logger = getLogger(__name__)

# interpolator types that implement fold constraints
FOLD_INTERPOLATOR_TYPES = ("DFI", "P2Fold")


class ModelFeatureFactory:
    @staticmethod
    def build_foliation(
        model: GeologicalModel,
        series_surface_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "FDI",
        nelements: int = LoopStructuralConfig.nelements,
        tol: float | None = None,
        faults: list[GeologicalFeature | str] | None = None,
        **kwargs: object,
    ) -> GeologicalFeature | None:
        """Build a foliation feature from model data.

        Parameters
        ----------
        series_surface_name : string
            corresponding to the feature_name in the data
        series_surface_data : pd.DataFrame, optional
            data frame containing the surface data
        interpolatortype : str
            the type of interpolator to use, default is 'FDI'
        nelements : int
            the number of elements to use in the series surface
        tol : float, optional
            tolerance for the solver, if not specified uses the model default
        faults : list, optional
            list of faults to be used in the series surface, if not specified uses the model faults
        kwargs

        Returns
        -------
        feature : GeologicalFeature
            the created geological feature

        Notes
        -----
        This function creates an instance of a
        :class:`LoopStructural.modelling.features.builders.GeologicalFeatureBuilder` and will return
        a :class:`LoopStructural.modelling.features.builders.GeologicalFeature`
        The feature is not interpolated until either
        :meth:`LoopStructural.modelling.features.builders.GeologicalFeature.evaluate_value` is called or
        :meth:`LoopStructural.modelling.core.GeologicalModel.update`

        An interpolator will be chosen by calling :meth:`LoopStructural.GeologicalModel.get_interpolator`

        """
        # if tol is not specified use the model default
        if tol is None:
            tol = model.tol

        series_builder = GeologicalFeatureBuilder(
            bounding_box=model.bounding_box,
            interpolatortype=interpolatortype,
            nelements=nelements,
            name=series_surface_name,
            model=model,
            **kwargs,
        )
        # add data
        if data is None:
            data = model.data.loc[model.data["feature_name"] == series_surface_name]

        if data.shape[0] == 0:
            logger.warning("No data for {series_surface_data}, skipping")
            return
        series_builder.add_data_from_data_frame(
            model.prepare_data(data, include_feature_name=False)
        )
        model._add_faults(series_builder, features=faults)

        # build feature
        series_feature = series_builder.feature
        series_builder.update_build_arguments(kwargs | {"domain": True, 'tol': tol})
        # this support is built for the entire model domain? Possibly would
        # could just pass a regular grid of points - mask by any above unconformities??

        series_feature.type = FeatureType.INTERPOLATED
        model._add_feature(series_feature, index=index)
        return series_feature

    @staticmethod
    def build_fold_frame(
        model: GeologicalModel,
        fold_frame_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "FDI",
        nelements: int = LoopStructuralConfig.nelements,
        tol: float | None = None,
        buffer: float = 0.1,
        **kwargs: object,
    ) -> FoldFrame | None:
        """Build a fold frame from model data.

        Parameters
        ----------
        fold_frame_name : string
            unique string in feature_name column
        fold_frame_data : pandas data frame
            if not specified uses the model data
        interpolatortype : str
            the type of interpolator to use, default is 'FDI'
        nelements : int
            the number of elements to use in the fold frame
        tol : float, optional
            tolerance for the solver
        buffer : float
            buffer to add to the bounding box of the fold frame
        **kwargs : dict
            additional parameters to be passed to the
            :class:`LoopStructural.modelling.features.builders.StructuralFrameBuilder`
            and :meth:`LoopStructural.modelling.features.builders.StructuralFrameBuilder.setup`
            and the interpolator, such as `domain` or `tol`


        Returns
        -------
        fold_frame : FoldFrame
            the created fold frame
        """
        if tol is None:
            tol = model.tol

        # create fault frame
        #
        fold_frame_builder = StructuralFrameBuilder(
            interpolatortype=interpolatortype,
            bounding_box=model.bounding_box.with_buffer(buffer),
            name=fold_frame_name,
            frame=FoldFrame,
            nelements=nelements,
            model=model,
            **kwargs,
        )
        # add data
        if data is None:
            data = model.data.loc[model.data["feature_name"] == fold_frame_name]
        if data.shape[0] == 0:
            logger.warning(f"No data for {fold_frame_name}, skipping")
            return
        fold_frame_builder.add_data_from_data_frame(
            model.prepare_data(data, include_feature_name=False)
        )
        model._add_faults(fold_frame_builder[0])
        model._add_faults(fold_frame_builder[1])
        model._add_faults(fold_frame_builder[2])
        kwargs["tol"] = tol
        fold_frame_builder.build(**kwargs)
        fold_frame = fold_frame_builder.frame

        fold_frame.type = FeatureType.STRUCTURALFRAME
        fold_frame.builder = fold_frame_builder
        model._add_feature(fold_frame, index=index)

        return fold_frame

    @staticmethod
    def build_folded_foliation(
        model: GeologicalModel,
        foliation_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "DFI",
        nelements: int = LoopStructuralConfig.nelements,
        buffer: float = 0.1,
        fold_frame: FoldFrame | None = None,
        svario: bool = True,
        tol: float | None = None,
        invert_fold_norm: bool = False,
        **kwargs: object,
    ) -> GeologicalFeature | None:
        """
        Create a folded foliation field from data and a fold frame.

        Parameters
        ----------
        foliation_data : str
            unique string in type column of data frame
        fold_frame :  FoldFrame
        svario  : Boolean
            whether to calculate svariograms, saves time if avoided
        interpolatortype : str
            ``"DFI"`` (default, piecewise linear) or ``"P2Fold"`` (piecewise quadratic).
            Other values are replaced by ``"DFI"`` with a warning.
        nelements : int
            requested structured-grid cell budget. Each cell splits into five
            tetrahedra. For ``"P2Fold"``, a request of 1e4 gave about 56k tets
            and 86k unknowns in the Bermagui validation; inspect the actual counts
        kwargs
            additional kwargs to be passed through to other functions

        Returns
        -------
        feature : GeologicalFeature
            created geological feature

        Notes
        -----
        - The P2 mesh is larger than the model box because of ``buffer``.
        - Building a folded foliation uses the fold interpolation code from Laurent et al., 2016
        and fold profile fitting from Grose et al., 2017. For more information about the fold modelling
        see :class:`LoopStructural.modelling.features.fold.FoldEvent`,
        :class:`LoopStructural.modelling.features.builders.FoldedFeatureBuilder`

        """
        if tol is None:
            tol = model.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = model.features[-1]
        if not isinstance(fold_frame, FoldFrame):
            raise TypeError("Please specify a FoldFrame")

        fold = FoldEvent(fold_frame, name=f"Fold_{foliation_name}", invert_norm=invert_fold_norm)

        if interpolatortype not in FOLD_INTERPOLATOR_TYPES:
            logger.warning(
                f"Folded foliation only supports {FOLD_INTERPOLATOR_TYPES} interpolators, "
                f"changing {interpolatortype} to DFI"
            )
            interpolatortype = "DFI"
        series_builder = FoldedFeatureBuilder(
            interpolatortype=interpolatortype,
            bounding_box=model.bounding_box.with_buffer(buffer),
            nelements=nelements,
            fold=fold,
            name=foliation_name,
            svario=svario,
            model=model,
            **kwargs,
        )
        if data is None:
            data = model.data.loc[model.data["feature_name"] == foliation_name]
        if data.shape[0] == 0:
            logger.warning(f"No data for {foliation_name}, skipping")
            return
        series_builder.add_data_from_data_frame(
            model.prepare_data(data, include_feature_name=False)
        )

        model._add_faults(series_builder)
        # build feature

        kwargs["tol"] = tol

        series_feature = series_builder.feature
        series_builder.update_build_arguments(kwargs)
        series_feature.type = FeatureType.FOLDED
        series_feature.fold = fold

        model._add_feature(series_feature, index)
        return series_feature

    @staticmethod
    def build_folded_fold_frame(
        model: GeologicalModel,
        fold_frame_name: str,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "FDI",
        nelements: int = LoopStructuralConfig.nelements,
        fold_frame: FoldFrame | None = None,
        tol: float | None = None,
        **kwargs: object,
    ) -> FoldFrame:
        """Build a folded fold frame from model data.

        Parameters
        ----------
        fold_frame_name : string
            name of the feature to be added
        fold_frame_data : pandas data frame, optional
            data frame containing the fold frame data, if not specified uses the model data
        interpolatortype : str
            ``"DFI"`` or ``"P2Fold"`` sets the interpolator of the first (fold) coordinate.
            Any other value gives ``"DFI"``. The other two coordinates always use FDI.
            ``nelements`` is a structured-grid cell budget, not a literal
            tetrahedron count. P2 adds edge nodes; inspect the actual dof.
        fold_frame : StructuralFrame, optional
            the fold frame for the fold if not specified uses last feature added
        nelements : int
            the number of elements to use in the fold frame
        tol : float, optional
            tolerance for the solver, if not specified uses the model default
        **kwargs : dict
            additional parameters to be passed to the
            :class:`LoopStructural.modelling.features.builders.StructuralFrameBuilder`
            and :meth:`LoopStructural.modelling.features.builders.StructuralFrameBuilder.setup`

        Returns
        -------
        fold_frame : FoldFrame
            created fold frame

        Notes
        -----
        This function build a structural frame where the first coordinate is constrained
        with a fold interpolator.
        Keyword arguments can be included to constrain

        - :meth:`LoopStructural.GeologicalModel.get_interpolator`
        - :class:`LoopStructural.StructuralFrameBuilder`
        - :meth:`LoopStructural.StructuralFrameBuilder.setup`
         - Building a folded foliation uses the fold interpolation code from Laurent et al., 2016
        and fold profile fitting from Grose et al., 2017. For more information about the fold modelling
        see :class:`LoopStructural.modelling.features.fold.FoldEvent`,
        :class:`LoopStructural.modelling.features.builders.FoldedFeatureBuilder`
        """
        if tol is None:
            tol = model.tol

        if fold_frame is None:
            logger.info("Using last feature as fold frame")
            fold_frame = model.features[-1]
        if not isinstance(fold_frame, FoldFrame):
            raise TypeError("Please specify a FoldFrame")
        fold = FoldEvent(fold_frame, name=f"Fold_{fold_frame_name}")

        # the first coordinate carries the fold constraints, the others are plain FDI
        interpolatortypes = [
            interpolatortype if interpolatortype in FOLD_INTERPOLATOR_TYPES else "DFI",
            "FDI",
            "FDI",
        ]
        fold_frame_builder = StructuralFrameBuilder(
            interpolatortype=interpolatortypes,
            bounding_box=model.bounding_box.with_buffer(kwargs.get("buffer", 0.1)),
            nelements=[nelements, nelements, nelements],
            name=fold_frame_name,
            fold=fold,
            frame=FoldFrame,
            model=model,
            **kwargs,
        )
        if data is None:
            data = model.data[model.data["feature_name"] == fold_frame_name]
        fold_frame_builder.add_data_from_data_frame(
            model.prepare_data(data, include_feature_name=False)
        )

        for i in range(3):
            model._add_faults(fold_frame_builder[i])
        # build feature
        kwargs["frame"] = FoldFrame
        kwargs["tol"] = tol
        fold_frame_builder.build(**kwargs)
        folded_fold_frame = fold_frame_builder.frame
        folded_fold_frame.builder = fold_frame_builder

        folded_fold_frame.type = FeatureType.STRUCTURALFRAME

        model._add_feature(folded_fold_frame, index=index)

        return folded_fold_frame

    @staticmethod
    def validate_intrusion_inputs(
        intrusion_name: str,
        intrusion_frame_name: str,
        intrusion_data: pd.DataFrame,
        intrusion_frame_data: pd.DataFrame,
        intrusion_frame_parameters: dict[str, object],
    ) -> None:
        """Fail fast at the intrusion creation boundary with a clear error.

        This names the missing piece instead of raising a bare `KeyError`
        several calls deep inside `IntrusionFrameBuilder`/`IntrusionBuilder`
        once building has already started. See ``INTRUSIONS.md`` finding 4.
        """
        if intrusion_data.empty:
            raise ValueError(
                f"No data found for intrusion '{intrusion_name}': check that "
                "model.data contains rows with feature_name == "
                f"'{intrusion_name}'"
            )
        if intrusion_frame_data.empty:
            raise ValueError(
                f"No data found for intrusion frame '{intrusion_frame_name}': "
                "check that model.data contains rows with feature_name == "
                f"'{intrusion_frame_name}'"
            )
        required_columns = ["intrusion_contact_type", "intrusion_side"]
        missing_columns = [c for c in required_columns if c not in intrusion_data.columns]
        if missing_columns:
            raise ValueError(
                f"Intrusion data for '{intrusion_name}' is missing required "
                f"column(s) {missing_columns}: 'intrusion_contact_type' marks "
                "each point as 'roof'/'top' or 'floor'/'base', and "
                "'intrusion_side' (boolean) marks points used to constrain "
                "the lateral extent"
            )
        contact_anisotropies = intrusion_frame_parameters.get("contact_anisotropies")
        if not contact_anisotropies:
            raise ValueError(
                "intrusion_frame_parameters['contact_anisotropies'] is "
                "required: provide a non-empty list of series-type features "
                "to use as the inflation-gradient proxy for the intrusion "
                "frame's coordinate 0"
            )

    @staticmethod
    def build_intrusion(
        model: GeologicalModel,
        intrusion_name: str,
        intrusion_frame_name: str,
        *,
        intrusion_frame_parameters: dict[str, object] | None = None,
        intrusion_lateral_extent_model: Callable[..., object] | None = None,
        intrusion_vertical_extent_model: Callable[..., object] | None = None,
        geometric_scaling_parameters: dict[str, object] | None = None,
        **kwargs: object,
    ) -> IntrusionFeature:
        """Build an intrusion feature from model data.

        Note
        -----
        An intrusion in built in two main steps:
        (1) Intrusion builder: intrusion builder creates the intrusion structural frame.
            This object is curvilinear coordinate system of the intrusion constrained with intrusion network points,
            and flow and inflation measurements (provided by the user).
            The intrusion network is a representation of the approximated location of roof or floor contact of the intrusion.
            This object might be constrained using the anisotropies of the host rock if the roof (or floor) contact is not well constrained.

        (2) Intrusion feature: simulation of lateral and vertical extent of intrusion within the model volume.
            The simulations outcome consist in thresholds distances along the structural frame coordinates
            that are used to constrained the extent of the intrusion.

        Parameters
        ----------
        intrusion_name :  string,
            name of intrusion feature in model data
        intrusion_frame_name :  string,
            name of intrusion frame in model data
        intrusion_lateral_extent_model = function,
            geometrical conceptual model for simulation of lateral extent
        intrusion_vertical_extent_model = function,
            geometrical conceptual model for simulation of vertical extent
        intrusion_frame_parameters = dictionary

        kwargs

        Returns
        -------
        intrusion feature

        """
        if intrusion_frame_parameters is None:
            intrusion_frame_parameters = {}
        if geometric_scaling_parameters is None:
            geometric_scaling_parameters = {}

        intrusion_data = model.data[model.data["feature_name"] == intrusion_name].copy()
        intrusion_frame_data = model.data[model.data["feature_name"] == intrusion_frame_name].copy()

        ModelFeatureFactory.validate_intrusion_inputs(
            intrusion_name,
            intrusion_frame_name,
            intrusion_data,
            intrusion_frame_data,
            intrusion_frame_parameters,
        )

        # -- get variables for intrusion frame interpolation
        gxxgz = kwargs.get("gxxgz", 0)
        gxxgy = kwargs.get("gxxgy", 0)
        gyxgz = kwargs.get("gyxgz", 0)

        interpolatortype = kwargs.get("interpolatortype", "PLI")
        # buffer = kwargs.get("buffer", 0.1)
        nelements = kwargs.get("nelements", LoopStructuralConfig.nelements)

        weights = [gxxgz, gxxgy, gyxgz]

        intrusion_frame_builder = IntrusionFrameBuilder(
            interpolatortype=interpolatortype,
            bounding_box=model.bounding_box.with_buffer(kwargs.get("buffer", 0.1)),
            nelements=kwargs.get("nelements", LoopStructuralConfig.nelements),
            name=intrusion_frame_name,
            model=model,
            **kwargs,
        )

        model._add_faults(intrusion_frame_builder)
        # intrusion_frame_builder.post_intrusion_faults = faults  # LG unused?

        # -- create intrusion frame using intrusion structures (steps and marginal faults) and flow/inflation measurements
        if len(intrusion_frame_parameters) == 0:
            logger.error("Please specify parameters to build intrusion frame")
        intrusion_frame_builder.set_intrusion_frame_parameters(
            intrusion_data, intrusion_frame_parameters
        )
        intrusion_frame_builder.create_constraints_for_c0()

        intrusion_frame_builder.set_intrusion_frame_data(intrusion_frame_data)

        ## -- create intrusion frame
        intrusion_frame_builder.build(
            nelements=nelements,
            w2=weights[0],
            w1=weights[1],
            gyxgz=weights[2],
        )

        intrusion_frame = intrusion_frame_builder.frame

        # -- create intrusion builder to compute distance thresholds along the frame coordinates
        intrusion_builder = IntrusionBuilder(
            intrusion_frame,
            model=model,
            # interpolator=interpolator,
            name=f"{intrusion_name}_feature",
            lateral_extent_model=intrusion_lateral_extent_model,
            vertical_extent_model=intrusion_vertical_extent_model,
            **kwargs,
        )
        intrusion_builder.set_data_for_extent_calculation(intrusion_data)

        intrusion_builder.update_build_arguments(
            {
                "geometric_scaling_parameters": geometric_scaling_parameters,
            }
        )

        intrusion_feature = intrusion_builder.feature
        model._add_feature(intrusion_feature)

        return intrusion_feature

    @staticmethod
    def build_domain_fault(
        model: GeologicalModel,
        fault_surface_data: pd.DataFrame | str,
        *,
        nelements: int = LoopStructuralConfig.nelements,
        interpolatortype: str = "FDI",
        index: int | None = None,
        **kwargs: object,
    ) -> UnconformityFeature:
        """Build a domain fault from model data.

        Parameters
        ----------
        fault_surface_data : string
            name of the domain fault data in the data frame

        Returns
        -------
        domain_Fault : GeologicalFeature
            the created domain fault

        Notes
        -----
        * :meth:`LoopStructural.GeologicalModel.get_interpolator`

        """
        domain_fault_feature_builder = GeologicalFeatureBuilder(
            bounding_box=model.bounding_box,
            interpolatortype=interpolatortype,
            nelements=nelements,
            name=fault_surface_data,
            model=model,
            **kwargs,
        )

        # add data
        unconformity_data = model.data.loc[model.data["feature_name"] == fault_surface_data]

        domain_fault_feature_builder.add_data_from_data_frame(unconformity_data)
        # look through existing features if there is a fault before an
        # unconformity
        # then add to the feature, once we get to an unconformity stop
        model._add_faults(domain_fault_feature_builder)

        # build feature
        domain_fault = domain_fault_feature_builder.feature
        domain_fault_feature_builder.update_build_arguments(kwargs)
        domain_fault.type = FeatureType.DOMAINFAULT
        model._add_feature(domain_fault, index=index)
        model._add_domain_fault_below(domain_fault)

        domain_fault_uc = UnconformityFeature(domain_fault, 0)
        # iterate over existing features and add the unconformity as a region
        # so the feature is only evaluated where the unconformity is positive
        return domain_fault_uc

    @staticmethod
    def build_fault(
        model: GeologicalModel,
        fault_name: str,
        displacement: float,
        *,
        index: int | None = None,
        data: pd.DataFrame | None = None,
        interpolatortype: str = "FDI",
        tol: float | None = None,
        fault_slip_vector: ArrayLike | None = None,
        fault_normal_vector: ArrayLike | None = None,
        fault_center: ArrayLike | None = None,
        major_axis: float | None = None,
        minor_axis: float | None = None,
        intermediate_axis: float | None = None,
        faultfunction: str = "BaseFault",
        faults: list[GeologicalFeature | str] | None = None,
        force_mesh_geometry: bool = False,
        points: bool = False,
        fault_buffer: float = 0.2,
        fault_trace_anisotropy: float = 0.0,
        fault_dip: float = 90,
        fault_dip_anisotropy: float = 0.0,
        fault_pitch: float | None = None,
        **kwargs: object,
    ) -> FaultSegment:
        """Build a fault from model data.

        Parameters
        ----------
        fault_name : string
            name of the fault surface data in the dataframe
        displacement : displacement magnitude
            displacement magnitude of the fault, in model units
        fault_data : pd.DataFrame, optional
            data frame containing the fault data, if not specified uses the model data
        major_axis : [type], optional
            [description], by default None
        minor_axis : [type], optional
            [description], by default None
        intermediate_axis : [type], optional
            [description], by default None
        kwargs : additional kwargs for Fault and interpolators

        Returns
        -------
        fault : FaultSegment
            created fault

        Notes
        -----
        * :meth:`LoopStructural.GeologicalModel.get_interpolator`
        * :class:`LoopStructural.modelling.features.builders.FaultBuilder`
        * :meth:`LoopStructural.modelling.features.builders.FaultBuilder.setup`
        """
        if faults is None:
            faults = []
        if "fault_extent" in kwargs and major_axis is None:
            major_axis = kwargs["fault_extent"]
        if "fault_influence" in kwargs and minor_axis is None:
            minor_axis = kwargs["fault_influence"]
        if "fault_vectical_radius" in kwargs and intermediate_axis is None:
            intermediate_axis = kwargs["fault_vectical_radius"]

        logger.info(f'Creating fault "{fault_name}"')
        logger.info(f"Displacement: {displacement}")
        logger.info(f"Tolerance: {tol}")
        logger.info(f"Fault function: {faultfunction}")
        logger.info(f"Fault slip vector: {fault_slip_vector}")
        logger.info(f"Fault center: {fault_center}")
        logger.info(f"Major axis: {major_axis}")
        logger.info(f"Minor axis: {minor_axis}")
        logger.info(f"Intermediate axis: {intermediate_axis}")
        if fault_slip_vector is not None:
            fault_slip_vector = np.array(fault_slip_vector, dtype="float")
        if fault_center is not None:
            fault_center = np.array(fault_center, dtype="float")

        for k, v in kwargs.items():
            logger.info(f"{k}: {v}")

        if tol is None:
            tol = model.tol
            # divide the tolerance by half of the minor axis, as this is the equivalent of the distance
            # of the unit vector
            # if minor_axis:
            # tol *= 0.1*minor_axis

        if displacement == 0:
            logger.warning(f"{fault_name} displacement is 0")

        if "data_region" in kwargs:
            kwargs.pop("data_region")
            logger.error("kwarg data_region currently not supported, disabling")
        displacement_scaled = displacement
        fault_frame_builder = FaultBuilder(
            interpolatortype,
            bounding_box=model.bounding_box,
            nelements=kwargs.pop("nelements", LoopStructuralConfig.nelements),
            name=fault_name,
            model=model,
            **kwargs,
        )
        if data is None:
            data = model.data.loc[model.data["feature_name"] == fault_name]
        if data.shape[0] == 0:
            logger.warning(f"No data for {fault_name}, skipping")
            return

        model._add_faults(fault_frame_builder, features=faults)
        # add data

        if fault_center is not None and ~np.isnan(fault_center).any():
            fault_center = model.scale(fault_center, inplace=False)
        # Keep the supplied fault-axis values unchanged; the previous self-assignment
        # was only present to satisfy a linter and did not affect behavior.
        fault_frame_builder.create_data_from_geometry(
            fault_frame_data=model.prepare_data(data, include_feature_name=False),
            fault_center=fault_center,
            fault_normal_vector=fault_normal_vector,
            fault_slip_vector=fault_slip_vector,
            minor_axis=minor_axis,
            major_axis=major_axis,
            intermediate_axis=intermediate_axis,
            points=points,
            force_mesh_geometry=force_mesh_geometry,
            fault_buffer=fault_buffer,
            fault_trace_anisotropy=fault_trace_anisotropy,
            fault_dip=fault_dip,
            fault_dip_anisotropy=fault_dip_anisotropy,
            fault_pitch=fault_pitch,
        )
        if "force_mesh_geometry" not in kwargs:
            fault_frame_builder.set_mesh_geometry(kwargs.get("fault_buffer", 0.2), 0)
        if "splay" in kwargs and "splayregion" in kwargs:
            fault_frame_builder.add_splay(kwargs["splay"], kwargs["splayregion"])

        kwargs["tol"] = tol
        fault_frame_builder.build(**kwargs)
        fault = fault_frame_builder.frame
        fault.displacement = displacement_scaled
        fault.faultfunction = faultfunction

        for f in reversed(model.features):
            if f.type == FeatureType.UNCONFORMITY:
                fault.add_region(f)
                break
        if displacement == 0:
            fault.type = FeatureType.INACTIVEFAULT
        model._add_feature(fault, index=index)

        return fault
