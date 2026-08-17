"""Recipe (JSON) and pickle serialization logic for GeologicalModel.

Extracted from GeologicalModel to separate serialization concerns from
feature-container orchestration (see API.md). GeologicalModel's
``@public_api``-decorated methods (``to_dict``, ``to_recipe_dict``,
``from_recipe_dict``, ``to_recipe_json``, ``from_recipe_json``,
``save_recipe``, ``load_recipe``, ``to_file``, ``from_file``) stay defined
directly on the class -- their ``__qualname__`` is part of the CI-checked
stable API surface (``tests/unit/test_public_api_contract.py``) -- and just
delegate to the staticmethods here.
"""

import json
import pathlib

import pandas as pd

from ...geometry import BoundingBox
from ...utils import LoopValueError, getLogger
from ..features import GeologicalFeature, StructuralFrame, UnconformityFeature
from ..features.fault import FaultSegment
from .stratigraphic_column import StratigraphicColumn

logger = getLogger(__name__)


class ModelSerializer:
    @staticmethod
    def feature_recipe_kind(feature):
        if isinstance(feature, GeologicalFeature):
            return "foliation"
        if isinstance(feature, StructuralFrame):
            return "structural_frame"
        if isinstance(feature, UnconformityFeature):
            return "unconformity"
        if isinstance(feature, FaultSegment):
            return "fault"
        return feature.__class__.__name__.lower()

    @staticmethod
    def to_dict(model):
        result = {}
        result["model"] = {}
        result["model"]["features"] = [f.name for f in model.features]
        result['model']['bounding_box'] = model.bounding_box.to_dict()
        result["model"]["stratigraphic_column"] = model.stratigraphic_column
        return result

    @staticmethod
    def to_recipe_dict(model, data_reference=None):
        recipe = {
            "schema": "LoopStructural.GeologicalModelRecipe",
            "version": 1,
            "model": {
                "bounding_box": model.bounding_box.to_dict(),
                "stratigraphic_column": model.stratigraphic_column.to_dict(),
                "data_source": None,
                "features": [],
            },
        }
        for feature in model.features:
            feature_entry = {
                "name": feature.name,
                "kind": ModelSerializer.feature_recipe_kind(feature),
                "faults": [
                    fault.name
                    for fault in getattr(feature, "faults", [])
                    if getattr(fault, "name", None)
                ],
                "regions": [],
            }
            recipe["model"]["features"].append(feature_entry)
        if data_reference is not None:
            recipe["model"]["data_source"] = {
                "kind": "reference",
                "path": str(pathlib.Path(data_reference)),
            }
        elif not model.data.empty:
            recipe["model"]["data_source"] = {
                "kind": "inline",
                "dataframe": model.data.to_dict(orient="split"),
            }
        return recipe

    @staticmethod
    def from_recipe_dict(cls, recipe):
        if not isinstance(recipe, dict):
            raise TypeError("recipe must be a dictionary")

        model_data = recipe.get("model", recipe)
        bounding_box = model_data.get("bounding_box")
        if isinstance(bounding_box, dict):
            bounding_box = BoundingBox.from_dict(bounding_box)
        if not isinstance(bounding_box, BoundingBox):
            raise TypeError("recipe must include a bounding_box dictionary")

        model = cls(bounding_box)

        data_source = model_data.get("data_source")
        if isinstance(data_source, dict):
            kind = data_source.get("kind")
            if kind == "reference":
                model.data = pd.read_csv(pathlib.Path(data_source["path"]))
            elif kind == "inline":
                dataframe = data_source.get("dataframe")
                if not isinstance(dataframe, dict):
                    raise TypeError("inline data_source must include a dataframe dictionary")
                model.data = pd.DataFrame(**dataframe)
            elif kind is not None:
                raise ValueError(f"Unsupported data_source kind: {kind}")
        elif isinstance(data_source, str):
            model.data = pd.read_csv(pathlib.Path(data_source))
        elif data_source is not None:
            raise TypeError("data_source must be a dictionary, string path, or None")

        stratigraphic_column = model_data.get("stratigraphic_column")
        if isinstance(stratigraphic_column, dict):
            model.stratigraphic_column = StratigraphicColumn.from_dict(stratigraphic_column)
        elif stratigraphic_column is not None:
            raise TypeError("stratigraphic_column must be a dictionary or None")

        features = model_data.get("features", [])
        if features is None:
            features = []
        if not isinstance(features, list):
            raise TypeError("features must be a list")

        feature_map = {}
        for feature_entry in features:
            if not isinstance(feature_entry, dict):
                raise TypeError("each feature entry must be a dictionary")
            feature_name = feature_entry.get("name")
            if not isinstance(feature_name, str):
                raise TypeError("each feature entry must include a string name")
            feature_data = model.data.loc[model.data["feature_name"] == feature_name].copy()
            if feature_data.empty:
                feature_data = None
            feature = model.create_and_add_foliation(feature_name, data=feature_data)
            if feature is None:
                raise ValueError(f"Could not recreate feature '{feature_name}' from recipe")
            feature_map[feature_name] = feature

        for feature_entry in features:
            feature_name = feature_entry.get("name")
            fault_names = feature_entry.get("faults", [])
            if not isinstance(fault_names, list):
                raise TypeError("faults must be a list")
            if fault_names:
                feature = feature_map[feature_name]
                feature.faults = [feature_map[name] for name in fault_names if name in feature_map]

        return model

    @staticmethod
    def to_recipe_json(model, data_reference=None, indent=2):
        recipe = ModelSerializer.to_recipe_dict(model, data_reference=data_reference)
        return json.dumps(recipe, indent=indent)

    @staticmethod
    def from_recipe_json(cls, json_str):
        if not isinstance(json_str, str):
            raise TypeError("json_str must be a string")
        try:
            recipe = json.loads(json_str)
        except json.JSONDecodeError as e:
            raise TypeError(f"json_str is not valid JSON: {e}")
        return ModelSerializer.from_recipe_dict(cls, recipe)

    @staticmethod
    def save_recipe(model, filename, data_reference=None):
        filename = pathlib.Path(filename)
        recipe = ModelSerializer.to_recipe_dict(model, data_reference=data_reference)
        with open(filename, "w") as f:
            json.dump(recipe, f, indent=2)
        logger.info(f"Recipe saved to {filename}")

    @staticmethod
    def load_recipe(cls, filename):
        filename = pathlib.Path(filename)
        if not filename.exists():
            raise FileNotFoundError(f"Recipe file not found: {filename}")
        with open(filename, "r") as f:
            recipe = json.load(f)
        logger.info(f"Recipe loaded from {filename}")
        return ModelSerializer.from_recipe_dict(cls, recipe)

    @staticmethod
    def to_file(model, file):
        try:
            import dill as pickle
        except ImportError:
            logger.error("Cannot write to file, dill not installed \n" "pip install dill")
            return
        try:
            logger.info(f"Writing GeologicalModel to: {file}")
            with open(file, "wb") as handle:
                pickle.dump(model, handle)
        except pickle.PicklingError:
            logger.error("Error saving file")

    @staticmethod
    def from_file(cls, file, allow_pickle: bool = True):
        if not allow_pickle:
            raise LoopValueError(
                "Pickle-based loading is disabled (allow_pickle=False). "
                f"Refusing to unpickle '{file}' because deserialising untrusted "
                "pickle/dill data can execute arbitrary code. If you generated "
                "this file yourself and trust its contents, call "
                "GeologicalModel.from_file(file, allow_pickle=True). Otherwise, "
                "use the JSON-based GeologicalModel.from_recipe_dict "
                "(paired with GeologicalModel.to_recipe_dict) as a safe "
                "alternative serialisation format."
            )
        logger.warning(
            f"Loading GeologicalModel from '{file}' using dill/pickle. "
            "Only load model files from trusted sources: deserialising a "
            "pickle file can execute arbitrary code. Pass allow_pickle=False "
            "to refuse pickle-based loading, or use "
            "GeologicalModel.from_recipe_dict for untrusted/JSON-based input."
        )
        try:
            import dill as pickle
        except ImportError:
            logger.error("Cannot import from file, dill not installed")
            return None
        path = pathlib.Path(file)
        if not path.is_file():
            raise LoopValueError(f"Cannot load model, file does not exist: {file}")
        try:
            with open(path, "rb") as f:
                model = pickle.load(f)
        except Exception as e:
            logger.error(f"Failed to load model from {file}: {e}")
            raise LoopValueError(f"Failed to load model from {file}: {e}") from e
        if isinstance(model, cls):
            logger.info("GeologicalModel initialised from file")
            return model
        else:
            logger.error(f"{file} does not contain a geological model")
            return None
