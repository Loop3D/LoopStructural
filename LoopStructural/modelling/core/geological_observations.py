"""
Geological Observations Data Structures

Provides intuitive, geologically-meaningful data structures for observations
that can be attached to geological objects in scenarios.
"""

from typing import Optional, List, Literal, Union
from dataclasses import dataclass, field
import numpy as np
import pandas as pd
from enum import Enum


class ObservationType(Enum):
    """Types of geological observations."""
    # Unit observations
    CONTACT = "contact"  # Points on unit boundary
    TANGENT = "tangent"  # Tangent vector to contact
    GRADIENT = "gradient"  # Gradient vector (normal to contact)
    ORIENTATION = "orientation"  # Strike/dip measurement
    INSIDE = "inside"  # Point known to be inside unit
    OUTSIDE = "outside"  # Point known to be outside unit
    ABOVE = "above"  # Point above the contact surface
    BELOW = "below"  # Point below the contact surface
    THICKNESS = "thickness"  # Thickness measurement
    
    # Fault observations
    TRACE = "trace"  # Fault trace point
    FAULT_ORIENTATION = "fault_orientation"  # Fault plane orientation
    DISPLACEMENT = "displacement"  # Displacement measurement
    HANGINGWALL = "hangingwall"  # Point in hanging wall
    FOOTWALL = "footwall"  # Point in footwall
    SLIP_VECTOR = "slip_vector"  # Slip direction vector
    
    # Fold observations
    AXIAL_SURFACE = "axial_surface"  # Point on axial surface
    FOLD_AXIS = "fold_axis"  # Fold axis orientation
    HINGE = "hinge"  # Point on fold hinge


@dataclass
class GeologicalObservation:
    """
    Base class for geological observations.
    
    All observations have a location and an observation type.
    Additional data depends on the specific observation.
    """
    location: np.ndarray  # [x, y, z]
    obs_type: ObservationType
    weight: float = 1.0
    comment: Optional[str] = None
    
    def __post_init__(self):
        self.location = np.asarray(self.location, dtype=float)
        if self.location.shape != (3,):
            raise ValueError("Location must be a 3D point [x, y, z]")


@dataclass
class ContactObservation(GeologicalObservation):
    """Point observation on a unit contact/boundary."""
    obs_type: ObservationType = field(default=ObservationType.CONTACT, init=False)
    scalar_value: Optional[float] = None  # Known scalar field value at this point


@dataclass
class OrientationObservation(GeologicalObservation):
    """
    Orientation observation (strike/dip or gradient vector).
    
    Can be specified as either:
    - Strike and dip (degrees)
    - Direct gradient vector
    - Tangent vector (will be converted to gradient)
    """
    obs_type: ObservationType = field(default=ObservationType.ORIENTATION, init=False)
    strike: Optional[float] = None  # Strike in degrees
    dip: Optional[float] = None  # Dip in degrees
    dip_direction: Optional[float] = None  # Dip direction in degrees
    gradient: Optional[np.ndarray] = None  # Direct gradient vector [gx, gy, gz]
    tangent: Optional[np.ndarray] = None  # Tangent vector (perpendicular to gradient)
    polarity: float = 1.0  # +1 or -1 for gradient direction
    
    def __post_init__(self):
        super().__post_init__()
        # Validate that we have either strike/dip or gradient
        has_strike_dip = self.strike is not None and self.dip is not None
        has_gradient = self.gradient is not None
        has_tangent = self.tangent is not None
        
        if not (has_strike_dip or has_gradient or has_tangent):
            raise ValueError("Must provide either strike/dip, gradient vector, or tangent vector")
        
        if has_gradient:
            self.gradient = np.asarray(self.gradient, dtype=float)
        if has_tangent:
            self.tangent = np.asarray(self.tangent, dtype=float)


@dataclass
class InsideOutsideObservation(GeologicalObservation):
    """
    Point known to be inside or outside a geological unit.
    Useful for constraining interpolation.
    """
    obs_type: ObservationType = field(default=ObservationType.INSIDE, init=False)
    inside: bool = True  # True if inside, False if outside


@dataclass
class AboveBelowObservation(GeologicalObservation):
    """
    Point known to be above or below a contact surface.
    """
    obs_type: ObservationType = field(default=ObservationType.ABOVE, init=False)
    above: bool = True  # True if above, False if below


@dataclass
class ThicknessObservation:
    """
    Thickness measurement for a unit.
    Not tied to a specific location but applies to the unit.
    """
    thickness: float
    location: Optional[np.ndarray] = None  # Optional location where measured
    weight: float = 1.0
    comment: Optional[str] = None


@dataclass
class FaultTraceObservation(GeologicalObservation):
    """Point on a fault trace (surface intersection)."""
    obs_type: ObservationType = field(default=ObservationType.TRACE, init=False)
    trace_direction: Optional[np.ndarray] = None  # Direction along trace


@dataclass
class FaultOrientationObservation(GeologicalObservation):
    """Fault plane orientation observation."""
    obs_type: ObservationType = field(default=ObservationType.FAULT_ORIENTATION, init=False)
    strike: Optional[float] = None
    dip: Optional[float] = None
    dip_direction: Optional[float] = None
    normal_vector: Optional[np.ndarray] = None  # Fault plane normal


@dataclass
class DisplacementObservation(GeologicalObservation):
    """Fault displacement measurement."""
    obs_type: ObservationType = field(default=ObservationType.DISPLACEMENT, init=False)
    displacement: Optional[float] = None  # Magnitude of displacement
    direction: Optional[np.ndarray] = None  # Displacement direction vector


@dataclass
class HangingwallFootwallObservation(GeologicalObservation):
    """Point in hanging wall or footwall of fault."""
    obs_type: ObservationType = field(default=ObservationType.HANGINGWALL, init=False)
    is_hangingwall: bool = True  # True for hanging wall, False for footwall


@dataclass
class SlipVectorObservation(GeologicalObservation):
    """Slip vector observation for fault."""
    obs_type: ObservationType = field(default=ObservationType.SLIP_VECTOR, init=False)
    slip_vector: Optional[np.ndarray] = None  # [sx, sy, sz]
    
    def __post_init__(self):
        super().__post_init__()
        if self.slip_vector is None:
            raise ValueError("slip_vector is required for SlipVectorObservation")
        self.slip_vector = np.asarray(self.slip_vector, dtype=float)


class ObservationCollection:
    """
    Collection of observations for a geological feature.
    
    Provides methods to add observations in a geologically intuitive way
    and convert them to LoopStructural's data format.
    """
    
    def __init__(self, feature_name: str):
        self.feature_name = feature_name
        self.observations: List[GeologicalObservation] = []
        self.thickness: Optional[float] = None
        
    def add_contact(self, location: Union[np.ndarray, List], 
                   weight: float = 1.0, scalar_value: Optional[float] = None,
                   comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a contact point observation."""
        obs = ContactObservation(
            location=location, 
            weight=weight,
            scalar_value=scalar_value,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_orientation(self, location: Union[np.ndarray, List],
                       strike: Optional[float] = None,
                       dip: Optional[float] = None,
                       dip_direction: Optional[float] = None,
                       gradient: Optional[np.ndarray] = None,
                       tangent: Optional[np.ndarray] = None,
                       polarity: float = 1.0,
                       weight: float = 1.0,
                       comment: Optional[str] = None) -> 'ObservationCollection':
        """Add an orientation observation (strike/dip or gradient)."""
        obs = OrientationObservation(
            location=location,
            strike=strike,
            dip=dip,
            dip_direction=dip_direction,
            gradient=gradient,
            tangent=tangent,
            polarity=polarity,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_inside_point(self, location: Union[np.ndarray, List],
                        weight: float = 1.0,
                        comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point known to be inside the unit."""
        obs = InsideOutsideObservation(
            location=location,
            inside=True,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_outside_point(self, location: Union[np.ndarray, List],
                         weight: float = 1.0,
                         comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point known to be outside the unit."""
        obs = InsideOutsideObservation(
            location=location,
            inside=False,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_above_point(self, location: Union[np.ndarray, List],
                       weight: float = 1.0,
                       comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point known to be above the contact surface."""
        obs = AboveBelowObservation(
            location=location,
            above=True,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_below_point(self, location: Union[np.ndarray, List],
                       weight: float = 1.0,
                       comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point known to be below the contact surface."""
        obs = AboveBelowObservation(
            location=location,
            above=False,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def set_thickness(self, thickness: float,
                     location: Optional[np.ndarray] = None) -> 'ObservationCollection':
        """Set the thickness for this unit."""
        self.thickness = thickness
        return self
    
    def add_fault_trace(self, location: Union[np.ndarray, List],
                       trace_direction: Optional[np.ndarray] = None,
                       weight: float = 1.0,
                       comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a fault trace point."""
        obs = FaultTraceObservation(
            location=location,
            trace_direction=trace_direction,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_fault_orientation(self, location: Union[np.ndarray, List],
                            strike: Optional[float] = None,
                            dip: Optional[float] = None,
                            dip_direction: Optional[float] = None,
                            normal_vector: Optional[np.ndarray] = None,
                            weight: float = 1.0,
                            comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a fault plane orientation."""
        obs = FaultOrientationObservation(
            location=location,
            strike=strike,
            dip=dip,
            dip_direction=dip_direction,
            normal_vector=normal_vector,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_displacement(self, location: Union[np.ndarray, List],
                        displacement: float,
                        direction: Optional[np.ndarray] = None,
                        weight: float = 1.0,
                        comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a displacement measurement."""
        obs = DisplacementObservation(
            location=location,
            displacement=displacement,
            direction=direction,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_hangingwall_point(self, location: Union[np.ndarray, List],
                             weight: float = 1.0,
                             comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point in the hanging wall."""
        obs = HangingwallFootwallObservation(
            location=location,
            is_hangingwall=True,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_footwall_point(self, location: Union[np.ndarray, List],
                          weight: float = 1.0,
                          comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a point in the footwall."""
        obs = HangingwallFootwallObservation(
            location=location,
            is_hangingwall=False,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def add_slip_vector(self, location: Union[np.ndarray, List],
                       slip_vector: np.ndarray,
                       weight: float = 1.0,
                       comment: Optional[str] = None) -> 'ObservationCollection':
        """Add a slip vector observation."""
        obs = SlipVectorObservation(
            location=location,
            slip_vector=slip_vector,
            weight=weight,
            comment=comment
        )
        self.observations.append(obs)
        return self
    
    def to_dataframe(self) -> pd.DataFrame:
        """
        Convert observations to LoopStructural DataFrame format.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with columns: feature_name, X, Y, Z, val, nx, ny, nz, 
            gx, gy, gz, tx, ty, tz, w, coord, polarity
        """
        rows = []
        
        for obs in self.observations:
            row = {
                'feature_name': self.feature_name,
                'X': obs.location[0],
                'Y': obs.location[1],
                'Z': obs.location[2],
                'val': np.nan,
                'nx': np.nan,
                'ny': np.nan,
                'nz': np.nan,
                'gx': np.nan,
                'gy': np.nan,
                'gz': np.nan,
                'tx': np.nan,
                'ty': np.nan,
                'tz': np.nan,
                'w': obs.weight,
                'coord': 0,
                'polarity': 1.0,
            }
            
            # Handle different observation types
            if isinstance(obs, ContactObservation):
                row['coord'] = 0  # Interface constraint
                if obs.scalar_value is not None:
                    row['val'] = obs.scalar_value
                    
            elif isinstance(obs, OrientationObservation):
                row['coord'] = 1  # Gradient constraint
                row['polarity'] = obs.polarity
                
                if obs.gradient is not None:
                    row['gx'] = obs.gradient[0]
                    row['gy'] = obs.gradient[1]
                    row['gz'] = obs.gradient[2]
                elif obs.tangent is not None:
                    row['tx'] = obs.tangent[0]
                    row['ty'] = obs.tangent[1]
                    row['tz'] = obs.tangent[2]
                # strike/dip will be converted by LoopStructural's prepare_data
                    
            elif isinstance(obs, InsideOutsideObservation):
                row['coord'] = 2  # Inequality constraint
                row['val'] = 1.0 if obs.inside else -1.0
                    
            elif isinstance(obs, AboveBelowObservation):
                row['coord'] = 2  # Inequality constraint
                row['val'] = 1.0 if obs.above else -1.0
                
            elif isinstance(obs, FaultTraceObservation):
                row['coord'] = 0  # On fault surface
                row['val'] = 0.0  # Fault scalar field = 0 on trace
                if obs.trace_direction is not None:
                    row['tx'] = obs.trace_direction[0]
                    row['ty'] = obs.trace_direction[1]
                    row['tz'] = obs.trace_direction[2]
                    
            elif isinstance(obs, FaultOrientationObservation):
                row['coord'] = 1  # Gradient/normal constraint
                if obs.normal_vector is not None:
                    row['gx'] = obs.normal_vector[0]
                    row['gy'] = obs.normal_vector[1]
                    row['gz'] = obs.normal_vector[2]
                # strike/dip handled by LoopStructural
                    
            elif isinstance(obs, DisplacementObservation):
                # Displacement observations are metadata, not point constraints
                # Store in val for now, will be handled specially
                row['val'] = obs.displacement
                if obs.direction is not None:
                    row['gx'] = obs.direction[0]
                    row['gy'] = obs.direction[1]
                    row['gz'] = obs.direction[2]
                    
            elif isinstance(obs, HangingwallFootwallObservation):
                row['coord'] = 2  # Inequality constraint
                row['val'] = 1.0 if obs.is_hangingwall else -1.0
                
            elif isinstance(obs, SlipVectorObservation):
                row['coord'] = 1  # Direction constraint
                row['gx'] = obs.slip_vector[0]
                row['gy'] = obs.slip_vector[1]
                row['gz'] = obs.slip_vector[2]
            
            rows.append(row)
        
        return pd.DataFrame(rows)
    
    def __len__(self):
        return len(self.observations)
    
    def __repr__(self):
        return f"ObservationCollection('{self.feature_name}', {len(self)} observations)"
