from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import Callable, Union

import numpy as np

DirectionProvider = Union[np.ndarray, Callable[[np.ndarray], np.ndarray]]


@dataclass(frozen=True)
class DirectionalRegularisation:
    weight: float
    direction: DirectionProvider
    name: str = "directional regularisation"


@dataclass(frozen=True)
class RegularisationConfig:
    isotropic: float | None = None
    directional: tuple[DirectionalRegularisation, ...] = ()


def _is_directional_mapping(value) -> bool:
    if not isinstance(value, dict):
        return False
    return any(key in value for key in ("direction", "vector", "weight", "name"))


def _coerce_directional_term(
    value, index: int = 0, default_name: str = "directional regularisation"
) -> DirectionalRegularisation:
    if isinstance(value, DirectionalRegularisation):
        return value

    if isinstance(value, dict):
        direction = value.get("direction", value.get("vector", None))
        if direction is None:
            raise ValueError("Directional regularisation entries require 'direction' or 'vector'")
        weight = value.get("weight", None)
        if weight is None:
            raise ValueError("Directional regularisation entries require 'weight'")
        name = value.get("name", f"{default_name} {index + 1}")
        return DirectionalRegularisation(
            weight=float(weight),
            direction=direction,
            name=name,
        )

    raise TypeError("Directional regularisation must be a DirectionalRegularisation or dict entry")


def coerce_directional_regularisation(
    value,
    default_name: str = "directional regularisation",
) -> tuple[DirectionalRegularisation, ...]:
    if value is None:
        return ()

    if isinstance(value, (DirectionalRegularisation, dict)):
        return (_coerce_directional_term(value, default_name=default_name),)

    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, np.ndarray)):
        if len(value) == 0:
            return ()
        if all(np.isscalar(item) for item in value):
            raise TypeError(
                "A sequence of scalars is not a valid directional regularisation config"
            )
        return tuple(
            _coerce_directional_term(item, index=index, default_name=default_name)
            for index, item in enumerate(value)
        )

    raise TypeError("Unsupported directional regularisation configuration")


def coerce_regularisation_config(
    regularisation=None,
    directional_regularisation=None,
) -> RegularisationConfig:
    isotropic = None
    directional = ()

    if isinstance(regularisation, RegularisationConfig):
        isotropic = regularisation.isotropic
        directional = regularisation.directional
    elif np.isscalar(regularisation) and regularisation is not None:
        isotropic = float(regularisation)
    elif _is_directional_mapping(regularisation):
        directional = coerce_directional_regularisation(regularisation)
    elif isinstance(regularisation, dict):
        isotropic_value = regularisation.get("isotropic", regularisation.get("weight", None))
        if (
            isotropic_value is not None
            and "direction" not in regularisation
            and "vector" not in regularisation
        ):
            isotropic = float(isotropic_value)
        directional = coerce_directional_regularisation(
            regularisation.get("directional", ()),
        )
    elif regularisation is not None:
        raise TypeError("Unsupported regularisation configuration")

    if directional_regularisation is not None:
        directional = directional + coerce_directional_regularisation(directional_regularisation)

    return RegularisationConfig(isotropic=isotropic, directional=directional)
