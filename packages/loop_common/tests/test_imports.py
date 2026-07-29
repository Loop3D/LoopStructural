import pytest


def test_import_common_modules():
    """Test if modules from the common package can be imported."""
    try:
        import loop_common.geometry
        import loop_common.io
        import loop_common.logging
        import loop_common.math
        import loop_common.supports
    except ImportError as e:
        pytest.fail(f"Failed to import a module from common: {e}")


def test_import_get_logger():
    """Test if get_logger can be imported from common.logging."""
    try:
        from loop_common.logging import get_logger

        assert callable(get_logger), "get_logger is not callable"
    except ImportError as e:
        pytest.fail(f"Failed to import get_logger from loop_common.logging: {e}")
