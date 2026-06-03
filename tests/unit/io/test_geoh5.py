import pytest
geoh5py = pytest.importorskip("geoh5py")
from LoopStructural.export.geoh5 import add_group_to_geoh5, add_points_to_geoh5, add_points_from_df

from pathlib import Path
from LoopStructural.datatypes import ValuePoints, VectorPoints
import numpy as np

@pytest.fixture
def tmp_path():
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        yield Path(tmpdir)
@pytest.fixture
def test_setup(tmp_path):
    filename = tmp_path / "test.geoh5"
    with geoh5py.workspace.Workspace.create(filename) as workspace:
        yield filename
        workspace.close()

def test_add_group_to_geoh5(test_setup):
    filename = test_setup
    group_uid = add_group_to_geoh5(filename, groupname="TestGroup")

    with geoh5py.workspace.Workspace(filename) as workspace:
        assert workspace.get_entity(group_uid)[0].name == "TestGroup"

def test_add_points_to_geoh5(test_setup):
    filename = test_setup
    group_uid = add_group_to_geoh5(filename, groupname="TestGroup")
    points = ValuePoints(
        name="TestPoints",
        locations=[[0, 0, 0], [1, 1, 1], [2, 2, 2]],
        values=[10., 20, 30],
    )
    add_points_to_geoh5(filename, points, groupname=group_uid)
    with geoh5py.workspace.Workspace(filename) as workspace:
        point_entity = workspace.get_entity("TestPoints")[0]
        assert point_entity.name == "TestPoints"
        assert point_entity.vertices.tolist() == [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        assert np.sum(point_entity.get_data("values")[0].values-np.array([10., 20., 30.])) == 0

def test_add_vector_points_to_geoh5(test_setup):
    filename = test_setup
    group_uid = add_group_to_geoh5(filename, groupname="TestGroup")
    points = VectorPoints(
        name="TestVectorPoints",
        locations=[[0, 0, 0], [1, 1, 1], [2, 2, 2]],
        vectors=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    )
    add_points_to_geoh5(filename, points, groupname=group_uid)
    with geoh5py.workspace.Workspace(filename) as workspace:
        point_entity = workspace.get_entity("TestVectorPoints")[0]
        assert point_entity.name == "TestVectorPoints"
        assert point_entity.vertices.tolist() == [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        assert np.sum(point_entity.get_data("vx")[0].values-np.array([1., 0., 0.])) == 0
        assert np.sum(point_entity.get_data("vy")[0].values-np.array([0., 1., 0.])) == 0
        assert np.sum(point_entity.get_data("vz")[0].values-np.array([0., 0., 1.])) == 0

def test_add_df_to_geoh5(test_setup):
    import pandas as pd
    filename = test_setup
    group_uid = add_group_to_geoh5(filename, groupname="TestGroup")
    df = pd.DataFrame({
        'X': [0, 1, 2],
        'Y': [0, 1, 2],
        'Z': [0, 1, 2],
        'value': [10., 20., 30.],
    })
    add_points_from_df(filename, df, name='df_points', groupname=group_uid)
    with geoh5py.workspace.Workspace(filename) as workspace:
        point_entity = workspace.get_entity("TestGroup")[0].children[0]
        assert point_entity.name == "df_points"
        assert point_entity.vertices.tolist() == [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        assert np.sum(point_entity.get_data("value")[0].values-np.array([10., 20., 30.])) == 0

def test_add_df_with_alternate_xyz_to_geoh5(test_setup):
    import pandas as pd
    filename = test_setup
    group_uid = add_group_to_geoh5(filename, groupname="TestGroup")
    df = pd.DataFrame({
        'EAST': [0, 1, 2],
        'NORTH': [0, 1, 2],
        'RL': [0, 1, 2],
        'value': [10., 20., 30.],
    })
    add_points_from_df(filename, df, name='df_points',groupname=group_uid, x_col='EAST', y_col='NORTH', z_col='RL')
    with geoh5py.workspace.Workspace(filename) as workspace:
        point_entity = workspace.get_entity("TestGroup")[0].children[0]
        assert point_entity.name == "df_points"
        assert point_entity.vertices.tolist() == [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        assert np.sum(point_entity.get_data("value")[0].values-np.array([10., 20., 30.])) == 0