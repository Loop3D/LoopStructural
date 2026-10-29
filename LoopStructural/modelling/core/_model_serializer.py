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

import numpy as np
import pandas as pd

from ...geometry import BoundingBox
from ...utils import LoopValueError, getLogger
from ..features import GeologicalFeature, StructuralFrame, UnconformityFeature
from ..features.fault import FaultSegment
from .stratigraphic_column import StratigraphicColumn

logger = getLogger(__name__)


# Parameters of create_and_add_feature that are not stored with a feature's creation record:
# ``index`` is the position (the recipe keeps the order), ``faults`` are other features (the recipe
# keeps their names), and ``data`` is a data frame (see ``custom_data``).
_NOT_RECORDED = ("index", "faults", "data")

# Feature types a recipe can rebuild from the parameters stored with it.
_SUPPORTED_FEATURE_TYPES = ("foliation",)


def _json_safe(value: object, where: str) -> object:
    """Return ``value`` as something ``json`` can write, or raise ``LoopValueError``."""
    if value is None or isinstance(value, (bool, int, float, str)):
        return value
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (list, tuple)):
        return [_json_safe(v, where) for v in value]
    if isinstance(value, dict):
        return {str(k): _json_safe(v, f"{where}[{k!r}]") for k, v in value.items()}
    raise LoopValueError(
        f"{where} is a {type(value).__name__}, which cannot be written to a recipe. "
        "Only numbers, strings, booleans, None, lists and dictionaries of these can be stored."
    )


class ModelSerializer:
    @staticmethod
    def feature_recipe_kind(feature):
        # subclasses first: FaultSegment is a StructuralFrame and UnconformityFeature a
        # GeologicalFeature
        if isinstance(feature, FaultSegment):
            return "fault"
        if isinstance(feature, UnconformityFeature):
            return "unconformity"
        if isinstance(feature, GeologicalFeature):
            return "foliation"
        if isinstance(feature, StructuralFrame):
            return "structural_frame"
        return feature.__class__.__name__.lower()

    @staticmethod
    def record_creation(feature: object, feature_type: str, params: dict) -> None:
        """Remember how ``feature`` was created, for the recipe.

        Called by ``GeologicalModel.create_and_add_feature`` with the parameters that were
        *requested* (for example ``nelements``, not the element count that was realised).
        """
        feature._recipe_creation = {
            "feature_type": feature_type,
            "parameters": {k: v for k, v in params.items() if k not in _NOT_RECORDED},
            "custom_data": params.get("data") is not None,
        }

    @staticmethod
    def feature_build_entry(feature: object) -> dict:
        """Return the ``build`` part of a feature's recipe entry.

        Raises
        ------
        LoopValueError
            if the feature cannot be rebuilt from a recipe yet
        """
        kind = ModelSerializer.feature_recipe_kind(feature)
        creation = getattr(feature, "_recipe_creation", None)
        if creation is None:
            raise LoopValueError(
                f"Feature '{feature.name}' ({kind}) cannot be stored in a recipe: there is no "
                "record of how it was created. Recipes store features created with "
                "GeologicalModel.create_and_add_feature (or a create_and_add_* method)."
            )
        feature_type = creation["feature_type"]
        if kind != "foliation" or feature_type not in _SUPPORTED_FEATURE_TYPES:
            raise LoopValueError(
                f"Feature '{feature.name}' ({kind}, created as '{feature_type}') cannot be "
                "stored in a recipe yet: recipes currently support foliations only."
            )
        if creation["custom_data"]:
            raise LoopValueError(
                f"Feature '{feature.name}' was created with its own data=, which a recipe cannot "
                "store: a recipe rebuilds a feature from the rows of the model data whose "
                "feature_name matches. Put the data in model.data instead."
            )
        parameters = {
            key: _json_safe(value, f"parameter '{key}' of feature '{feature.name}'")
            for key, value in creation["parameters"].items()
        }
        return {"feature_type": feature_type, "parameters": parameters}

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
                "build": ModelSerializer.feature_build_entry(feature),
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
            build = feature_entry.get("build")
            if build is None:
                # a recipe written before build parameters were stored
                kind = feature_entry.get("kind", "foliation")
                if kind != "foliation":
                    raise LoopValueError(
                        f"Cannot rebuild feature '{feature_name}' of kind '{kind}' from this "
                        "recipe: recipes currently support foliations only."
                    )
                logger.warning(
                    f"The recipe has no build parameters for feature '{feature_name}' (it was "
                    "written by an older version); building it with the default parameters"
                )
                build = {"feature_type": "foliation", "parameters": {}}
            if not isinstance(build, dict) or not isinstance(build.get("parameters", {}), dict):
                raise TypeError("a feature's build entry must be a dictionary with parameters")
            feature = model.create_and_add_feature(
                build.get("feature_type", "foliation"),
                feature_name,
                **build.get("parameters", {}),
            )
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
