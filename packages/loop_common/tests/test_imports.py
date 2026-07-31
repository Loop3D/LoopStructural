import pytest


def test_import_common_modules():
    """Test if modules from the common package can be imported."""
    try:
        import loop_common.geometry
        import loop_common.io
        import loop_common.logging
        import loop_common.math  # noqa: F401
    except ImportError as e:
        pytest.fail(f"Failed to import a module from common: {e}")


def test_import_get_logger():
    """Test if get_logger can be imported from common.logging."""
    try:
        from loop_common.logging import get_logger

        assert callable(get_logger), "get_logger is not callable"
    except ImportError as e:
        pytest.fail(f"Failed to import get_logger from loop_common.logging: {e}")


def test_import_shared_geometry_and_utils():
    """Test if geometry and utils helpers are exposed from loop_common."""
    try:
        from loop_common.geometry import BoundingBox, Surface
        from loop_common.utils import getLogger, rng

        assert callable(getLogger)
        assert hasattr(rng, "random")
        assert BoundingBox is not None
        assert Surface is not None
    except ImportError as e:
        pytest.fail(f"Failed to import shared geometry/utils from loop_common: {e}")
