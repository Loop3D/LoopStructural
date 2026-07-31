"""Guards the "stable" API surface documented in API.md.

Signatures of every @public_api(tier="stable")-decorated method are
snapshotted in tests/fixtures/api_surface_snapshot.json. This test fails if
that surface changes (added, removed, or a signature edited) without the
change being reflected in both the snapshot and COMPAT.md, per the
versioning policy in ROADMAP.md.
"""

import json
from pathlib import Path

import LoopStructural.modelling.core.geological_model  # noqa: F401  (registers @public_api entries)
from LoopStructural.utils import get_stable_surface

SNAPSHOT_PATH = Path(__file__).parents[1] / "fixtures" / "api_surface_snapshot.json"
COMPAT_PATH = Path(__file__).parents[2] / "COMPAT.md"


def _load_snapshot():
    return json.loads(SNAPSHOT_PATH.read_text())


def test_stable_surface_matches_snapshot_or_is_logged_in_compat():
    snapshot = _load_snapshot()
    current = get_stable_surface()
    compat_text = COMPAT_PATH.read_text() if COMPAT_PATH.exists() else ""

    added = sorted(set(current) - set(snapshot))
    removed = sorted(set(snapshot) - set(current))
    changed = sorted(
        name
        for name in set(current) & set(snapshot)
        if current[name] != snapshot[name]
    )

    undocumented = [
        name
        for name in removed + changed
        if name.split(".")[-1] not in compat_text
    ]

    assert not undocumented, (
        "Stable API surface changed without a COMPAT.md entry for: "
        f"{undocumented}. Removed: {removed}. Changed: {changed}."
    )
    assert not added, (
        "New stable API methods are not yet captured in "
        f"{SNAPSHOT_PATH}: {added}. Add them to the snapshot once the "
        "signature is considered final."
    )


def test_stable_surface_is_non_empty():
    assert len(get_stable_surface()) > 0
