import uuid
import numpy as np
from datetime import datetime
from typing import Annotated, Any, Optional
from pydantic import BaseModel, Field, ConfigDict, PlainSerializer, BeforeValidator, TypeAdapter
from pathlib import Path
from loop_common.logging import get_logger as getLogger

logger = getLogger(__name__)
# --- 1. The NumPy Type Logic ---


def validate_numpy(v: Any) -> np.ndarray:
    """Ensures input is converted to a numpy array."""
    if isinstance(v, np.ndarray):
        return v
    try:
        return np.array(v)
    except Exception as e:
        raise ValueError(f"Could not convert {type(v)} to numpy array") from e


# Define a 'NumpyArray' type that:
# - Converts lists/tuples to arrays during input (BeforeValidator)
# - Converts arrays to lists during JSON export (PlainSerializer)
NumpyArray = Annotated[
    np.ndarray,
    BeforeValidator(validate_numpy),
    PlainSerializer(lambda x: x.tolist(), return_type=list),
]

# --- 2. The Base Entity ---


class LoopEntity(BaseModel):
    """
    The atomic building block for all Loop objects.
    Provides identity, validation, and serialization.
    """

    # Allow Pydantic to handle non-pydantic types (like numpy arrays)
    model_config = ConfigDict(
        arbitrary_types_allowed=True,
        validate_assignment=True,  # Validate if user changes a value later
        extra="forbid",  # Prevent accidental typos from creating new fields
    )

    uuid: str = Field(
        default_factory=lambda: str(uuid.uuid4()), description="Permanent unique identifier"
    )

    name: Optional[str] = Field(default=None, description="Human-readable label")

    last_modified: str = Field(
        default_factory=lambda: datetime.now().isoformat(),
        description="ISO timestamp of last change",
    )

    def mark_modified(self):
        """Manually trigger a timestamp update."""
        self.last_modified = datetime.now().isoformat()

    @classmethod
    def from_json(cls, json_str: str):
        """Helper to reconstruct the object from a JSON string."""
        return cls.model_validate_json(json_str)

    def to_json(self, indent: int = 2) -> str:
        """Helper to export to JSON string."""
        return self.model_dump_json(indent=indent)

    def to_yaml(self) -> str:
        """Export to YAML string."""
        try:
            import yaml
        except ImportError as exc:
            raise ImportError("PyYAML is required for YAML export: pip install pyyaml") from exc
        return yaml.dump(self.model_dump(mode="json"), sort_keys=False, allow_unicode=True)

    def to_dict(self) -> dict:
        """Export to a plain dict (JSON-safe: numpy arrays converted to lists)."""
        return self.model_dump(mode="json")

    def save(self, filename: Path):
        filename = Path(filename)
        filetype = filename.suffix.lstrip(".")
        if filetype == "json":
            with open(filename,'w') as f:
                f.write(self.to_json())
        elif filetype == 'yaml':
            with open(filename,'w') as f:
                f.write(self.to_yaml())
        else:
            logger.warning(f"Unknown filetype {filetype}. Not saving.")
