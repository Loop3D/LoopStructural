"""Extension point backing ``GeologicalModel.create_and_add_feature`` (see ``API.md``).

Maps a feature-type string to a factory callable ``factory(model, name,
**params) -> feature``. New feature types (e.g. a future intrusion-workflow
rewrite) register a factory here instead of requiring changes to
``GeologicalModel``'s source.
"""

from typing import Callable, Dict, List


class FeatureBuilderRegistry:
    _factories: Dict[str, Callable] = {}

    @classmethod
    def register(cls, feature_type: str, factory: Callable) -> None:
        cls._factories[feature_type] = factory

    @classmethod
    def create(cls, feature_type: str, model, name: str, **params):
        if feature_type not in cls._factories:
            raise ValueError(
                f"Unknown feature_type '{feature_type}'. Registered types: "
                f"{cls.registered_types()}"
            )
        return cls._factories[feature_type](model, name, **params)

    @classmethod
    def registered_types(cls) -> List[str]:
        return sorted(cls._factories)
