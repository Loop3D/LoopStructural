import itertools
from math import factorial

import numpy as np
import pytest

from loop_interpolation._quadrature import tet_quadrature, tet_quadrature_points

RULES = ["centroid", "gauss4"]


def _exact_mean(a):
    # mean of r^a0 s^a1 t^a2 over the reference tet: 3! * prod(a_i!) / (3 + sum a)!
    return 6.0 * np.prod([factorial(i) for i in a]) / factorial(3 + sum(a))


@pytest.mark.parametrize("rule", RULES)
def test_weights_sum_to_one(rule):
    bary, w = tet_quadrature(rule)
    assert np.isclose(w.sum(), 1.0)
    assert np.allclose(bary.sum(axis=1), 1.0)


@pytest.mark.parametrize(
    "rule,max_degree", [("centroid", 1), ("gauss4", 2)]
)
def test_exact_for_polynomials(rule, max_degree):
    bary, w = tet_quadrature(rule)
    rst = bary[:, 1:]
    for a in itertools.product(range(max_degree + 1), repeat=3):
        if sum(a) > max_degree:
            continue
        f = np.prod(rst ** np.array(a), axis=1)
        assert np.isclose(w @ f, _exact_mean(a)), a


def test_unknown_rule():
    with pytest.raises(ValueError):
        tet_quadrature("nope")
    with pytest.raises(ValueError):
        tet_quadrature_points(np.zeros((4, 3)), np.arange(4)[None], "nope")


def test_mapping_indices_and_points():
    nodes = np.array(
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 1]], dtype=float
    )
    elements = np.array([[0, 1, 2, 3], [1, 2, 3, 4]])
    pts, idx, w = tet_quadrature_points(nodes, elements, "gauss4")
    assert pts.shape == (8, 3)
    assert np.array_equal(idx, [0, 0, 0, 0, 1, 1, 1, 1])
    assert np.allclose(w.reshape(2, 4).sum(axis=1), 1.0)
    for e in range(2):
        mean = (w[idx == e, None] * pts[idx == e]).sum(axis=0)
        assert np.allclose(mean, nodes[elements[e]].mean(axis=0))
    pc, ic, wc = tet_quadrature_points(nodes, elements, "centroid")
    assert np.array_equal(ic, [0, 1])
    assert np.allclose(pc, nodes[elements].mean(axis=1))
