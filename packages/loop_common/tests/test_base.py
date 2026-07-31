import time
import uuid
from datetime import datetime

import numpy as np
import pytest
import yaml
from loop_common.base import LoopEntity, NumpyArray
from pydantic import ConfigDict, TypeAdapter, ValidationError

# --- construction / defaults ---


def test_default_uuid_is_valid_and_unique():
    e1 = LoopEntity()
    e2 = LoopEntity()
    uuid.UUID(e1.uuid)
    uuid.UUID(e2.uuid)
    assert e1.uuid != e2.uuid


def test_default_name_is_none():
    e = LoopEntity()
    assert e.name is None


def test_explicit_name_is_kept():
    e = LoopEntity(name="test")
    assert e.name == "test"


def test_default_last_modified_is_recent_iso_timestamp():
    before = datetime.now()
    e = LoopEntity()
    after = datetime.now()
    ts = datetime.fromisoformat(e.last_modified)
    assert before <= ts <= after


def test_explicit_uuid_is_kept():
    fixed = str(uuid.uuid4())
    e = LoopEntity(uuid=fixed)
    assert e.uuid == fixed


# --- mark_modified ---


def test_mark_modified_updates_timestamp_only():
    e = LoopEntity(name="test")
    old_ts = datetime.fromisoformat(e.last_modified)
    old_uuid = e.uuid
    old_name = e.name

    time.sleep(0.01)
    e.mark_modified()

    new_ts = datetime.fromisoformat(e.last_modified)
    assert new_ts > old_ts
    assert e.uuid == old_uuid
    assert e.name == old_name


# --- validation behaviour (extra="forbid", validate_assignment) ---


def test_unknown_field_is_rejected():
    with pytest.raises(ValidationError):
        LoopEntity(name="test", unknown_field=1)


def test_assignment_is_validated():
    e = LoopEntity(name="test")
    e.name = "renamed"
    assert e.name == "renamed"

    with pytest.raises(ValidationError):
        e.name = 123


# --- serialization: json ---


def test_json_roundtrip():
    e = LoopEntity(name="roundtrip")
    j = e.to_json()
    e2 = LoopEntity.from_json(j)
    assert e2.uuid == e.uuid
    assert e2.name == e.name
    assert e2.last_modified == e.last_modified


def test_to_json_is_indentable():
    e = LoopEntity(name="test")
    compact = e.to_json(indent=None)
    indented = e.to_json(indent=2)
    assert "\n" not in compact
    assert "\n" in indented


# --- serialization: dict / yaml ---


def test_to_dict_is_json_safe():
    d = LoopEntity(name="test").to_dict()
    assert d["name"] == "test"
    assert isinstance(d["uuid"], str)
    assert isinstance(d["last_modified"], str)


def test_to_yaml_roundtrips_via_yaml_load():
    e = LoopEntity(name="yaml-test")
    loaded = yaml.safe_load(e.to_yaml())
    assert loaded["uuid"] == e.uuid
    assert loaded["name"] == "yaml-test"
    assert loaded["last_modified"] == e.last_modified


# --- save() ---


def test_save_json_writes_loadable_file(tmp_path):
    e = LoopEntity(name="saved")
    target = tmp_path / "entity.json"
    e.save(target)

    assert target.exists()
    loaded = LoopEntity.from_json(target.read_text())
    assert loaded.uuid == e.uuid
    assert loaded.name == "saved"


def test_save_yaml_writes_loadable_file(tmp_path):
    e = LoopEntity(name="saved-yaml")
    target = tmp_path / "entity.yaml"
    e.save(target)

    assert target.exists()
    loaded = yaml.safe_load(target.read_text())
    assert loaded["uuid"] == e.uuid
    assert loaded["name"] == "saved-yaml"


def test_save_unknown_filetype_does_not_write(tmp_path):
    e = LoopEntity(name="unsaved")
    target = tmp_path / "entity.txt"
    e.save(target)
    assert not target.exists()


# --- NumpyArray type ---


def test_numpyarray_typeadapter_validate_and_dump():
    ta = TypeAdapter(NumpyArray, config=ConfigDict(arbitrary_types_allowed=True))
    arr = ta.validate_python([1, 2, 3])
    assert isinstance(arr, np.ndarray)
    dumped = ta.dump_python(arr)
    assert dumped == [1, 2, 3]


def test_numpyarray_passthrough_for_existing_array():
    ta = TypeAdapter(NumpyArray, config=ConfigDict(arbitrary_types_allowed=True))
    original = np.array([1.0, 2.0, 3.0])
    validated = ta.validate_python(original)
    assert validated is original


def test_numpyarray_rejects_unconvertible_input():
    ta = TypeAdapter(NumpyArray, config=ConfigDict(arbitrary_types_allowed=True))
    with pytest.raises(ValidationError):
        ta.validate_python([[1, 2], [3, 4, 5]])
