"""Compatibility re-export: the vector maths helpers now live in loop_common."""

from loop_common.math._maths import (
    azimuthplunge2vector,
    dipdipdirection2vector,
    get_dip_vector,
    get_strike_vector,
    get_vectors,
    gradient_from_tetrahedron,
    normal_vector_to_dip_and_dip_direction,
    normal_vector_to_strike_and_dip,
    plungeazimuth2vector,
    regular_tetraherdron_for_points,
    rotate,
    rotation,
    strikedip2vector,
)

__all__ = [
    "azimuthplunge2vector",
    "dipdipdirection2vector",
    "get_dip_vector",
    "get_strike_vector",
    "get_vectors",
    "gradient_from_tetrahedron",
    "normal_vector_to_dip_and_dip_direction",
    "normal_vector_to_strike_and_dip",
    "plungeazimuth2vector",
    "regular_tetraherdron_for_points",
    "rotate",
    "rotation",
    "strikedip2vector",
]
