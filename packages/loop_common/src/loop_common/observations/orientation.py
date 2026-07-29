import numpy as np
from loop_common.base import LoopEntity, NumpyArray
from loop_common.math import strikedip2vector, dipdipdirection2vector, plungeazimuth2vector
from pydantic import model_validator
from enum import Enum


class OrientationType(str, Enum):
    PLANE = "plane"
    LINEATION = "lineation"
    TANGENT = "tangent"


class OrientationObservation(LoopEntity):
    """Strike/Dip or Dip/DipDirection measurements."""

    coords: NumpyArray
    vector: NumpyArray  # Normal vector to the surface
    magnitude: NumpyArray
    polarity: NumpyArray  # 1 for upright, -1 for overturned
    type: OrientationType

    @model_validator(mode="after")
    def check_dimensions(self):
        if self.coords.shape[-1] != 3:
            raise ValueError("Coords must have shape (N, 3)")
        if self.vector.shape[-1] != 3:
            raise ValueError("Vector must have shape (N, 3)")
        if self.magnitude.shape[0] != self.coords.shape[0]:
            raise ValueError("Magnitude must have same length as coords")
        if self.polarity.shape[0] != self.coords.shape[0]:
            raise ValueError("Polarity must have same length as coords")
        if self.vector.shape[0] != self.coords.shape[0]:
            raise ValueError("Vector must have same length as coords")
        return self

    @classmethod
    def from_strike_dip(
        cls,
        coords: np.ndarray,
        strike: np.ndarray,
        dip: np.ndarray,
        polarity: np.ndarray,
        name: str | None = None,
    ):
        """Create an OrientationObservation from strike/dip measurements."""
        # Convert strike/dip to normal vector
        # This is a simplified conversion assuming right-hand rule and that strike is measured clockwise from north
        vector = strikedip2vector(strike, dip)
        magnitude = np.ones_like(strike)  # Placeholder for magnitude, could be set to

        return cls(
            name=name,
            coords=coords,
            vector=vector,
            magnitude=magnitude,
            polarity=polarity,
            type=OrientationType.PLANE,
        )

    @classmethod
    def from_dip_direction_and_dip(
        cls,
        coords: np.ndarray,
        dip_direction: np.ndarray,
        dip: np.ndarray,
        polarity: np.ndarray | None = None,
        name: str | None = None,
    ):
        """Create an OrientationObservation from dip direction/dip measurements."""
        # Convert dip direction/dip to normal vector
        if not hasattr(dip_direction, 'len'):
            dip_direction = np.ones_like(coords[:, 0]) * dip_direction
        if not hasattr(dip, 'len'):
            dip = np.ones_like(coords[:, 0]) * dip
        vector = dipdipdirection2vector(dip_direction, dip)
        magnitude = np.ones_like(
            dip_direction
        )  # Placeholder for magnitude, could be set to something else
        if polarity is None:
            polarity = np.ones_like(dip_direction)  # Default to upright if not provided
        return cls(
            name=name,
            coords=coords,
            vector=vector,
            magnitude=magnitude,
            polarity=polarity,
            type=OrientationType.PLANE,
        )

    @classmethod
    def from_plunge_and_plunge_direction(
        cls,
        coords: np.ndarray,
        plunge_direction: np.ndarray,
        plunge: np.ndarray,
        polarity: np.ndarray,
        name: str | None = None,
    ):
        """Create an OrientationObservation from plunge direction/plunge measurements."""
        # Convert plunge direction/plunge to normal vector
        vector = plungeazimuth2vector(plunge, plunge_direction)
        magnitude = np.ones_like(
            plunge_direction
        )  # Placeholder for magnitude, could be set to something else

        return cls(
            name=name,
            coords=coords,
            vector=vector,
            magnitude=magnitude,
            polarity=polarity,
            type=OrientationType.PLANE,
        )


# Backwards-compatible alias expected by other modules
Orientation = OrientationObservation
