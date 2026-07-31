import numpy as np
import pytest

from LoopStructural.geometry import BoundingBox, Surface
from LoopStructural.utils._surface import LoopIsosurfacer


@pytest.fixture()
def small_bbox():
    return BoundingBox(origin=[0, 0, 0], maximum=[1, 1, 1], nsteps=[10, 10, 10])


def plane_field(xyz):
    """Scalar field equal to (x - 0.5): zero isosurface is the x=0.5 plane."""
    xyz = np.asarray(xyz)
    return xyz[:, 0] - 0.5


class MockInterpolator:
    def evaluate_value(self, xyz):
        return plane_field(xyz)


def test_requires_interpolator_or_callable(small_bbox):
    with pytest.raises(ValueError):
        LoopIsosurfacer(small_bbox)


def test_cannot_specify_both_interpolator_and_callable(small_bbox):
    with pytest.raises(ValueError):
        LoopIsosurfacer(small_bbox, interpolator=MockInterpolator(), callable=plane_field)


def test_constructor_uses_interpolator_evaluate_value(small_bbox):
    iso = LoopIsosurfacer(small_bbox, interpolator=MockInterpolator())
    assert callable(iso.callable)
    # calling it should behave the same as calling evaluate_value directly
    xyz = small_bbox.regular_grid()
    assert np.allclose(iso.callable(xyz), plane_field(xyz))


def test_fit_extracts_isosurface_at_given_value(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=[0.0], name="plane")

    assert len(surfaces) == 1
    surface = surfaces[0]
    assert isinstance(surface, Surface)
    assert surface.name == "plane"
    # all vertices should lie approximately on the x=0.5 plane
    assert np.allclose(surface.vertices[:, 0], 0.5, atol=1e-6)
    assert np.allclose(surface.values, 0.0)
    assert surface.triangles.shape[1] == 3


def test_fit_with_single_value_and_list_name_uses_individual_name(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=[0.0], name=["custom_name"])

    assert len(surfaces) == 1
    assert surfaces[0].name == "custom_name"


def test_fit_with_multiple_values_names_include_isovalue(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=[-0.25, 0.25], name="iso")

    assert len(surfaces) == 2
    names = sorted(s.name for s in surfaces)
    assert names == sorted(["iso_-0.25", "iso_0.25"])


def test_fit_with_none_uses_mean_value(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=None)

    assert len(surfaces) == 1
    # mean of min/max of (x - 0.5) over the bounding box grid is ~0
    assert np.allclose(surfaces[0].vertices[:, 0], 0.5, atol=1e-6)


def test_fit_with_int_generates_multiple_evenly_spaced_isosurfaces(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=3, name="multi")

    assert len(surfaces) == 3
    x_values = sorted(s.vertices[0, 0] for s in surfaces)
    # evenly spaced with a 5% buffer inside [-0.5, 0.5] range of plane_field
    assert x_values[0] < x_values[1] < x_values[2]


def test_fit_with_int_less_than_one_raises(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    with pytest.raises(ValueError):
        iso.fit(values=-1)


def test_fit_assigns_colours(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    surfaces = iso.fit(values=[-0.25, 0.25], name="iso", colours=["red", "blue"])

    colours = {s.colour for s in surfaces}
    assert colours == {"red", "blue"}


def test_fit_skips_isovalue_outside_of_data_range(small_bbox):
    iso = LoopIsosurfacer(small_bbox, callable=plane_field)
    # plane_field ranges roughly from -0.5 to 0.5 over the bounding box,
    # so a value well outside that range cannot be marched and should be
    # skipped (with a warning) rather than raising.
    surfaces = iso.fit(values=[100.0])
    assert surfaces == []
