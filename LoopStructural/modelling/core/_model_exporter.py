"""Surface/block-model export logic for GeologicalModel (see API.md).

Extracted from GeologicalModel to separate export/visualization data prep
from feature-container orchestration. GeologicalModel's @public_api methods
(``get_fault_surfaces``, ``get_stratigraphic_surfaces``, ``get_block_model``,
``save``) stay defined directly on the class -- their __qualname__ is part
of the CI-checked stable API surface -- and delegate to the staticmethods
here.
"""

from __future__ import annotations

import pathlib
from typing import TYPE_CHECKING

from ...geometry import StructuredGrid
from ...utils import getLogger, surface_list

if TYPE_CHECKING:
    from .geological_model import GeologicalModel

logger = getLogger(__name__)


class ModelExporter:
    @staticmethod
    def get_fault_surfaces(model: GeologicalModel, faults: list[str] | None = None) -> surface_list:
        if faults is None:
            faults = []
        surfaces = []
        if len(faults) == 0:
            faults = model.fault_names()

        for f in faults:
            surfaces.extend(model.get_feature_by_name(f).surfaces([0], model.bounding_box))
        return surfaces

    @staticmethod
    def get_stratigraphic_surfaces(
        model: GeologicalModel, units: list[str] | None = None, bottoms: bool = True
    ) -> surface_list:
        if units is None:
            units = []
        ## TODO change the stratigraphic column to its own class and have methods to get the relevant surfaces
        surfaces = []
        units = []
        if model.stratigraphic_column is None:
            return []
        units = model.stratigraphic_column.get_isovalues()
        units_for_group = {}
        for name, u in units.items():
            if u['group'] not in model:
                logger.warning(f"Group {u['group']} not found in model")
                continue
            if u['group'] not in units_for_group:
                units_for_group[u['group']] = []
            u['name'] = name
            units_for_group[u['group']].append(u)
        for group, us in units_for_group.items():
            feature = model.get_feature_by_name(group)
            values = [u['value'] for u in us]
            colours = [u['colour'] for u in us]
            names = [u['name'] for u in us]
            surfaces.extend(
                feature.surfaces(values, model.bounding_box, name=names, colours=colours)
            )

        return surfaces

    @staticmethod
    def get_block_model(
        model: GeologicalModel, name: str = 'block model'
    ) -> tuple[StructuredGrid, list]:
        # NOTE: bounding_box.structured_grid() returns loop_common's
        # interpolation-support StructuredGrid (no properties dict); use
        # LoopStructural's own geometry StructuredGrid for storing values.
        grid = StructuredGrid(
            origin=model.bounding_box.origin,
            step_vector=model.bounding_box.step_vector,
            nsteps=model.bounding_box.nsteps,
            name=name,
        )

        grid.cell_properties['stratigraphy'] = model.evaluate_model(
            model.rescale(model.bounding_box.cell_centres())
        )
        return grid, model.stratigraphic_ids()

    @staticmethod
    def save(
        model: GeologicalModel,
        filename: str,
        block_model: bool = True,
        stratigraphic_surfaces: bool = True,
        fault_surfaces: bool = True,
        stratigraphic_data: bool = True,
        fault_data: bool = True,
    ) -> None:
        path = pathlib.Path(filename)
        extension = path.suffix
        parent = path.parent
        name = path.stem
        stratigraphic_surfaces = model.get_stratigraphic_surfaces()
        if fault_surfaces:
            for s in model.get_fault_surfaces():
                ## geoh5 can save everything into the same file
                if extension == ".geoh5" or extension == '.omf':
                    s.save(filename)
                else:
                    s.save(f'{parent}/{name}_{s.name}{extension}')
        if stratigraphic_surfaces:
            for s in model.get_stratigraphic_surfaces():
                if extension == ".geoh5" or extension == '.omf':
                    s.save(filename)
                else:
                    s.save(f'{parent}/{name}_{s.name}{extension}')
        if block_model:
            grid, _ids = model.get_block_model()
            if extension == ".geoh5" or extension == '.omf':
                grid.save(filename)
            else:
                grid.save(f'{parent}/{name}_block_model{extension}')
        if stratigraphic_data and model.stratigraphic_column is not None:
            for group in model.stratigraphic_column:
                if group == "faults":
                    continue
                for data in model.__getitem__(group).get_data():
                    if extension == ".geoh5" or extension == '.omf':
                        data.save(filename)
                    else:
                        data.save(f'{parent}/{name}_{group}_data{extension}')
        if fault_data:
            for f in model.fault_names():
                for d in model.__getitem__(f).get_data():
                    if extension == ".geoh5" or extension == '.omf':

                        d.save(filename)
                    else:
                        d.save(f'{parent}/{name}_{group}{extension}')
