"""Guards the parts of the stable surface (API.md) that the signature
snapshot (test_public_api_contract.py) can't check on its own: that
documented module paths stay importable, that documented top-level symbols
still exist, and that documented Enum members aren't silently renamed or
removed (their names/values are part of the contract too -- e.g. persisted
in `GeologicalModel.to_dict()`/recipe JSON, or compared directly by callers).

`PLUGIN_MODULE_PATHS` mirrors the "QGIS-plugin compatibility" list in
ROADMAP.md verbatim -- keep the two in sync. This runs locally with plain
`pytest`, unlike the fuller check in `.github/workflows/qgis-compat.yml`,
which additionally runs the plugin's own test suite against this branch but
needs that repo checked out.
"""

import importlib

import pytest

PLUGIN_MODULE_PATHS = [
    "LoopStructural.modelling.core.fault_topology",
    "LoopStructural.modelling.features",
    "LoopStructural.modelling.features.fold",
    "LoopStructural.modelling.features.builders",
    "LoopStructural.modelling.features._feature_converters",
    "LoopStructural.modelling.core.stratigraphic_column",
    "LoopStructural.datatypes",
    "LoopStructural.utils",
]


@pytest.mark.parametrize("module_path", PLUGIN_MODULE_PATHS)
def test_plugin_relied_on_module_path_importable(module_path):
    importlib.import_module(module_path)


def test_plugin_relied_on_top_level_symbols_importable():
    from LoopStructural import (  # noqa: F401
        GeologicalModel,
        FaultTopology,
        StratigraphicColumn,
        getLogger,
    )


def test_documented_stable_classes_importable():
    """Classes API.md lists as stable regardless of GeologicalModel usage."""
    from LoopStructural.modelling.core.fault_topology import (  # noqa: F401
        FaultRelationshipType,
    )
    from LoopStructural.modelling.core.stratigraphic_column import (  # noqa: F401
        StratigraphicColumnElementType,
    )
    from LoopStructural.modelling.features import (  # noqa: F401
        FeatureType,
        StructuralFrame,
    )
    from LoopStructural.modelling.features.fold import FoldFrame  # noqa: F401
    from LoopStructural.modelling.features.builders import (  # noqa: F401
        StructuralFrameBuilder,
        FaultBuilder,
        GeologicalFeatureBuilder,
        FoldedFeatureBuilder,
    )
    from LoopStructural.geometry import (  # noqa: F401
        BoundingBox,
        Surface,
        ValuePoints,
        VectorPoints,
    )
    from LoopStructural.utils.observer import Observable  # noqa: F401


# Enum members as of the "Stable surface" section in API.md (2026-07-30).
# Renaming/removing a member is breaking (values may be persisted in
# `to_dict()`/recipe JSON, or compared directly: `x.type == FeatureType.FAULT`).
# Adding new members is not breaking, so this only checks a subset.
PROTECTED_ENUM_MEMBERS = {
    "LoopStructural.modelling.core.fault_topology.FaultRelationshipType": [
        "ABUTTING",
        "FAULTED",
        "NONE",
    ],
    "LoopStructural.modelling.core.stratigraphic_column.StratigraphicColumnElementType": [
        "UNIT",
        "UNCONFORMITY",
    ],
    "LoopStructural.modelling.features.FeatureType": [
        "BASE",
        "INTERPOLATED",
        "STRUCTURALFRAME",
        "REGION",
        "FOLDED",
        "ANALYTICAL",
        "LAMBDA",
        "UNCONFORMITY",
        "INTRUSION",
        "FAULT",
        "DOMAINFAULT",
        "INACTIVEFAULT",
        "ONLAPUNCONFORMITY",
    ],
}


@pytest.mark.parametrize("qualname", PROTECTED_ENUM_MEMBERS)
def test_enum_members_not_removed(qualname):
    module_path, enum_name = qualname.rsplit(".", 1)
    enum_cls = getattr(importlib.import_module(module_path), enum_name)
    current_members = {member.name for member in enum_cls}
    missing = set(PROTECTED_ENUM_MEMBERS[qualname]) - current_members
    assert not missing, f"{qualname} is missing members: {sorted(missing)}"
