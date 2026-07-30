"""Registry backing the API contract documented in ``API.md``.

``@public_api`` is a no-op at call time; it only records the decorated
callable's qualified name, signature, and stability tier so that
``tests/unit/test_public_api_contract.py`` can snapshot the "stable" tier
and fail CI if it drifts without a matching ``COMPAT.md`` entry.
"""

import functools
import inspect
from typing import Callable, Dict, Literal

Tier = Literal["stable", "provisional"]

_REGISTRY: Dict[str, Dict[str, str]] = {}


def public_api(tier: Tier = "stable") -> Callable:
    def decorator(func: Callable) -> Callable:
        _REGISTRY[func.__qualname__] = {
            "tier": tier,
            "signature": str(inspect.signature(func)),
        }

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            return func(*args, **kwargs)

        return wrapper

    return decorator


def register_external_stable(qualname: str, obj: Callable, tier: Tier = "stable") -> None:
    """Register a class/function LoopStructural re-exports but doesn't define.

    `@public_api` can't be applied at the definition site for symbols owned by
    a separate package (e.g. `loop_common`'s `BoundingBox`/`Observable`) --
    that package has its own release cycle and shouldn't import LoopStructural
    internals. Call this instead, from the LoopStructural module that
    re-exports the symbol, to capture the same (qualname, signature, tier)
    entry for the snapshot test.
    """
    _REGISTRY[qualname] = {
        "tier": tier,
        "signature": str(inspect.signature(obj)),
    }


def get_registry() -> Dict[str, Dict[str, str]]:
    return dict(_REGISTRY)


def get_stable_surface() -> Dict[str, str]:
    return {
        name: entry["signature"]
        for name, entry in _REGISTRY.items()
        if entry["tier"] == "stable"
    }
