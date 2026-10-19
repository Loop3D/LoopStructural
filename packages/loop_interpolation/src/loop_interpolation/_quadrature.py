"""Tetrahedron volume quadrature rules.

Barycentric convention: c0 = 1 - r - s - t, c1..c3 = r, s, t.
Weights are normalised to sum to 1, so a rule gives the mean of a
function over the tetrahedron. Multiply by the element volume for
the integral.
"""

import numpy as np

_GAUSS4_A = 0.1381966011250105
_GAUSS4_B = 0.5854101966249685

_RULES = ("centroid", "gauss4")


def tet_quadrature(rule):
    """Return a tetrahedron quadrature rule.

    Parameters
    ----------
    rule : str
        'centroid' (1 point, exact for linear functions) or
        'gauss4' (4 points, exact for quadratic functions).

    Returns
    -------
    bary : np.ndarray
        (n, 4) barycentric coordinates of the points.
    weights : np.ndarray
        (n,) weights that sum to 1.

    Raises
    ------
    ValueError
        If the rule name is not known.
    """
    if rule == "centroid":
        return np.full((1, 4), 0.25), np.ones(1)
    if rule == "gauss4":
        a, b = _GAUSS4_A, _GAUSS4_B
        bary = np.full((4, 4), a)
        bary[np.arange(4), np.arange(4)] = b
        return bary, np.full(4, 0.25)
    raise ValueError(f"Unknown quadrature rule '{rule}', use one of {_RULES}")


def tet_quadrature_points(nodes, elements, rule):
    """Map a quadrature rule onto tetrahedra.

    Parameters
    ----------
    nodes : np.ndarray
        (n_nodes, 3) node coordinates.
    elements : np.ndarray
        (m, >=4) element node indices. Only the 4 corners are used.
    rule : str
        Rule name, see `tet_quadrature`.

    Returns
    -------
    points : np.ndarray
        (m * n, 3) points, grouped by element.
    tet_idx : np.ndarray
        (m * n,) index of the element (row of `elements`) for each point.
    weights : np.ndarray
        (m * n,) rule weights repeated for each element (sum to 1 per element).
    """
    bary, w = tet_quadrature(rule)
    corners = np.asarray(nodes)[np.asarray(elements)[:, :4]]  # (m, 4, 3)
    points = np.einsum("pk,mkd->mpd", bary, corners).reshape(-1, 3)
    n = len(w)
    tet_idx = np.repeat(np.arange(corners.shape[0]), n)
    weights = np.tile(w, corners.shape[0])
    return points, tet_idx, weights
