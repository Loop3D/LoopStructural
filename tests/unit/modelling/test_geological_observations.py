"""
Tests for geological observations data structures.
"""

import pytest
import numpy as np
import pandas as pd
from LoopStructural.modelling.core.geological_observations import (
    ObservationCollection,
    ContactObservation,
    OrientationObservation,
    InsideOutsideObservation,
    AboveBelowObservation,
    FaultTraceObservation,
    FaultOrientationObservation,
    DisplacementObservation,
    HangingwallFootwallObservation,
    SlipVectorObservation,
    ObservationType
)


class TestObservationDataClasses:
    """Test individual observation data classes."""
    
    def test_contact_observation(self):
        """Test ContactObservation creation."""
        obs = ContactObservation(location=[100, 200, 50], weight=1.5, scalar_value=10.0)
        assert obs.obs_type == ObservationType.CONTACT
        assert np.allclose(obs.location, [100, 200, 50])
        assert obs.weight == 1.5
        assert obs.scalar_value == 10.0
    
    def test_orientation_with_strike_dip(self):
        """Test OrientationObservation with strike/dip."""
        obs = OrientationObservation(
            location=[100, 200, 50],
            strike=45,
            dip=30
        )
        assert obs.obs_type == ObservationType.ORIENTATION
        assert obs.strike == 45
        assert obs.dip == 30
    
    def test_orientation_with_gradient(self):
        """Test OrientationObservation with gradient vector."""
        gradient = np.array([0.1, 0.2, 0.9])
        obs = OrientationObservation(
            location=[100, 200, 50],
            gradient=gradient
        )
        assert obs.obs_type == ObservationType.ORIENTATION
        assert np.allclose(obs.gradient, gradient)
    
    def test_orientation_requires_input(self):
        """Test that OrientationObservation requires either strike/dip or gradient."""
        with pytest.raises(ValueError, match="Must provide either"):
            OrientationObservation(location=[100, 200, 50])
    
    def test_inside_outside_observation(self):
        """Test InsideOutsideObservation."""
        inside_obs = InsideOutsideObservation(location=[100, 200, 50], inside=True)
        assert inside_obs.inside is True
        
        outside_obs = InsideOutsideObservation(location=[100, 200, 50], inside=False)
        assert outside_obs.inside is False
    
    def test_above_below_observation(self):
        """Test AboveBelowObservation."""
        above_obs = AboveBelowObservation(location=[100, 200, 50], above=True)
        assert above_obs.above is True
        
        below_obs = AboveBelowObservation(location=[100, 200, 50], above=False)
        assert below_obs.above is False
    
    def test_fault_trace_observation(self):
        """Test FaultTraceObservation."""
        obs = FaultTraceObservation(
            location=[100, 200, 0],
            trace_direction=np.array([1, 0, 0])
        )
        assert obs.obs_type == ObservationType.TRACE
        assert np.allclose(obs.trace_direction, [1, 0, 0])
    
    def test_fault_orientation_observation(self):
        """Test FaultOrientationObservation."""
        obs = FaultOrientationObservation(
            location=[100, 200, 50],
            strike=90,
            dip=60
        )
        assert obs.obs_type == ObservationType.FAULT_ORIENTATION
        assert obs.strike == 90
        assert obs.dip == 60
    
    def test_displacement_observation(self):
        """Test DisplacementObservation."""
        obs = DisplacementObservation(
            location=[100, 200, 50],
            displacement=100.0,
            direction=np.array([0, 0, -1])
        )
        assert obs.obs_type == ObservationType.DISPLACEMENT
        assert obs.displacement == 100.0
        assert np.allclose(obs.direction, [0, 0, -1])
    
    def test_hangingwall_footwall_observation(self):
        """Test HangingwallFootwallObservation."""
        hw_obs = HangingwallFootwallObservation(
            location=[100, 200, 50],
            is_hangingwall=True
        )
        assert hw_obs.is_hangingwall is True
        
        fw_obs = HangingwallFootwallObservation(
            location=[100, 200, 50],
            is_hangingwall=False
        )
        assert fw_obs.is_hangingwall is False
    
    def test_slip_vector_observation(self):
        """Test SlipVectorObservation."""
        slip = np.array([0, 0, -1])
        obs = SlipVectorObservation(
            location=[100, 200, 50],
            slip_vector=slip
        )
        assert obs.obs_type == ObservationType.SLIP_VECTOR
        assert np.allclose(obs.slip_vector, slip)
    
    def test_slip_vector_requires_vector(self):
        """Test that SlipVectorObservation requires slip_vector."""
        with pytest.raises(ValueError, match="slip_vector is required"):
            SlipVectorObservation(location=[100, 200, 50], slip_vector=None)


class TestObservationCollection:
    """Test ObservationCollection class."""
    
    def test_create_collection(self):
        """Test creating an empty collection."""
        collection = ObservationCollection("unit1")
        assert collection.feature_name == "unit1"
        assert len(collection) == 0
    
    def test_add_contact(self):
        """Test adding contact observations."""
        collection = ObservationCollection("unit1")
        collection.add_contact([100, 200, 50], weight=1.5, comment="Drillhole")
        
        assert len(collection) == 1
        assert isinstance(collection.observations[0], ContactObservation)
        assert collection.observations[0].weight == 1.5
    
    def test_add_orientation(self):
        """Test adding orientation observations."""
        collection = ObservationCollection("unit1")
        collection.add_orientation([100, 200, 50], strike=45, dip=30)
        
        assert len(collection) == 1
        assert isinstance(collection.observations[0], OrientationObservation)
    
    def test_add_inside_outside_points(self):
        """Test adding inside/outside observations."""
        collection = ObservationCollection("unit1")
        collection.add_inside_point([100, 200, 50])
        collection.add_outside_point([200, 300, 60])
        
        assert len(collection) == 2
        assert collection.observations[0].inside is True
        assert collection.observations[1].inside is False
    
    def test_add_above_below_points(self):
        """Test adding above/below observations."""
        collection = ObservationCollection("unit1")
        collection.add_above_point([100, 200, 60])
        collection.add_below_point([100, 200, 40])
        
        assert len(collection) == 2
        assert collection.observations[0].above is True
        assert collection.observations[1].above is False
    
    def test_add_fault_observations(self):
        """Test adding fault-specific observations."""
        collection = ObservationCollection("fault1")
        
        collection.add_fault_trace([100, 200, 0])
        collection.add_fault_orientation([100, 200, 50], strike=90, dip=60)
        collection.add_hangingwall_point([110, 200, 50])
        collection.add_footwall_point([90, 200, 50])
        collection.add_slip_vector([100, 200, 50], slip_vector=[0, 0, -1])
        
        assert len(collection) == 5
        assert isinstance(collection.observations[0], FaultTraceObservation)
        assert isinstance(collection.observations[1], FaultOrientationObservation)
        assert isinstance(collection.observations[2], HangingwallFootwallObservation)
        assert isinstance(collection.observations[3], HangingwallFootwallObservation)
        assert isinstance(collection.observations[4], SlipVectorObservation)
    
    def test_fluent_interface(self):
        """Test method chaining (fluent interface)."""
        collection = ObservationCollection("unit1")
        result = (collection
                  .add_contact([100, 200, 50])
                  .add_orientation([100, 200, 50], strike=45, dip=30)
                  .add_inside_point([150, 250, 50])
                  .set_thickness(25.0))
        
        assert result is collection  # Returns self
        assert len(collection) == 3
        assert collection.thickness == 25.0
    
    def test_to_dataframe(self):
        """Test conversion to DataFrame."""
        collection = ObservationCollection("unit1")
        collection.add_contact([100, 200, 50])
        collection.add_orientation([100, 200, 50], strike=45, dip=30)
        
        df = collection.to_dataframe()
        
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 2
        assert 'feature_name' in df.columns
        assert 'X' in df.columns
        assert 'Y' in df.columns
        assert 'Z' in df.columns
        assert all(df['feature_name'] == 'unit1')
        assert df.iloc[0]['coord'] == 0  # Contact
        assert df.iloc[1]['coord'] == 1  # Orientation
    
    def test_contact_to_dataframe(self):
        """Test contact observation DataFrame conversion."""
        collection = ObservationCollection("unit1")
        collection.add_contact([100, 200, 50], scalar_value=10.0, weight=2.0)
        
        df = collection.to_dataframe()
        
        assert df.iloc[0]['X'] == 100
        assert df.iloc[0]['Y'] == 200
        assert df.iloc[0]['Z'] == 50
        assert df.iloc[0]['val'] == 10.0
        assert df.iloc[0]['w'] == 2.0
        assert df.iloc[0]['coord'] == 0
    
    def test_orientation_to_dataframe(self):
        """Test orientation observation DataFrame conversion."""
        collection = ObservationCollection("unit1")
        gradient = [0.1, 0.2, 0.9]
        collection.add_orientation([100, 200, 50], gradient=gradient, polarity=-1.0)
        
        df = collection.to_dataframe()
        
        assert df.iloc[0]['coord'] == 1
        assert df.iloc[0]['polarity'] == -1.0
        assert df.iloc[0]['gx'] == 0.1
        assert df.iloc[0]['gy'] == 0.2
        assert df.iloc[0]['gz'] == 0.9
    
    def test_inequality_to_dataframe(self):
        """Test inequality constraint DataFrame conversion."""
        collection = ObservationCollection("unit1")
        collection.add_inside_point([100, 200, 50])
        collection.add_outside_point([200, 300, 60])
        collection.add_above_point([100, 200, 70])
        collection.add_below_point([100, 200, 30])
        
        df = collection.to_dataframe()
        
        assert all(df['coord'] == 2)  # All inequality constraints
        assert df.iloc[0]['val'] == 1.0  # Inside
        assert df.iloc[1]['val'] == -1.0  # Outside
        assert df.iloc[2]['val'] == 1.0  # Above
        assert df.iloc[3]['val'] == -1.0  # Below
    
    def test_fault_trace_to_dataframe(self):
        """Test fault trace DataFrame conversion."""
        collection = ObservationCollection("fault1")
        collection.add_fault_trace([100, 200, 0], trace_direction=[1, 0, 0])
        
        df = collection.to_dataframe()
        
        assert df.iloc[0]['coord'] == 0  # On surface
        assert df.iloc[0]['val'] == 0.0  # Fault surface
        assert df.iloc[0]['tx'] == 1.0
        assert df.iloc[0]['ty'] == 0.0
        assert df.iloc[0]['tz'] == 0.0
    
    def test_fault_orientation_to_dataframe(self):
        """Test fault orientation DataFrame conversion."""
        collection = ObservationCollection("fault1")
        normal = [1, 0, 0]
        collection.add_fault_orientation([100, 200, 50], normal_vector=normal)
        
        df = collection.to_dataframe()
        
        assert df.iloc[0]['coord'] == 1  # Gradient constraint
        assert df.iloc[0]['gx'] == 1.0
        assert df.iloc[0]['gy'] == 0.0
        assert df.iloc[0]['gz'] == 0.0
    
    def test_slip_vector_to_dataframe(self):
        """Test slip vector DataFrame conversion."""
        collection = ObservationCollection("fault1")
        slip = [0, 0, -1]
        collection.add_slip_vector([100, 200, 50], slip_vector=slip)
        
        df = collection.to_dataframe()
        
        assert df.iloc[0]['coord'] == 1  # Direction constraint
        assert df.iloc[0]['gx'] == 0.0
        assert df.iloc[0]['gy'] == 0.0
        assert df.iloc[0]['gz'] == -1.0
    
    def test_repr(self):
        """Test string representation."""
        collection = ObservationCollection("unit1")
        collection.add_contact([100, 200, 50])
        collection.add_contact([200, 300, 60])
        
        repr_str = repr(collection)
        assert "unit1" in repr_str
        assert "2 observations" in repr_str


class TestObservationIntegration:
    """Test integration with geological scenario."""
    
    def test_multiple_features(self):
        """Test observations for multiple features."""
        unit1 = ObservationCollection("unit1")
        unit1.add_contact([100, 200, 50])
        unit1.add_orientation([100, 200, 50], strike=45, dip=30)
        
        unit2 = ObservationCollection("unit2")
        unit2.add_contact([100, 200, 75])
        unit2.add_orientation([100, 200, 75], strike=50, dip=25)
        
        fault1 = ObservationCollection("fault1")
        fault1.add_fault_trace([100, 200, 0])
        fault1.add_fault_orientation([100, 200, 50], strike=90, dip=60)
        
        # Combine all dataframes
        df1 = unit1.to_dataframe()
        df2 = unit2.to_dataframe()
        df3 = fault1.to_dataframe()
        
        combined = pd.concat([df1, df2, df3], ignore_index=True)
        
        assert len(combined) == 6
        assert set(combined['feature_name'].unique()) == {'unit1', 'unit2', 'fault1'}
