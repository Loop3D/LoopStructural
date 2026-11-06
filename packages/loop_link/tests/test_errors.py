from __future__ import annotations

import dataclasses

import pytest

from loop_link.errors import (
    BusyError,
    InvalidInputError,
    LoopLinkError,
    MissingExtraError,
    MissingFeatureError,
    NoViewerError,
    ObjectTooLargeError,
    Refusal,
    TokenError,
    ViewerClosedError,
    ViewerRefusedError,
)

SIMPLE_ERRORS = [
    NoViewerError,
    TokenError,
    BusyError,
    InvalidInputError,
    ViewerClosedError,
]


@pytest.mark.parametrize("error_class", SIMPLE_ERRORS)
def test_simple_error_is_loop_link_error(error_class: type[LoopLinkError]) -> None:
    assert issubclass(error_class, LoopLinkError)
    assert issubclass(LoopLinkError, Exception)


@pytest.mark.parametrize("error_class", SIMPLE_ERRORS)
def test_simple_error_default_text_is_not_empty(error_class: type[LoopLinkError]) -> None:
    assert str(error_class()) != ""


@pytest.mark.parametrize("error_class", SIMPLE_ERRORS)
def test_simple_error_keeps_custom_text(error_class: type[LoopLinkError]) -> None:
    assert str(error_class("custom text")) == "custom text"


def test_missing_feature_error_keeps_feature_and_names_it() -> None:
    error = MissingFeatureError("viewer")

    assert isinstance(error, LoopLinkError)
    assert error.feature == "viewer"
    assert "viewer" in str(error)


def test_object_too_large_error_keeps_attributes_and_text() -> None:
    error = ObjectTooLargeError("mesh-1", 2048, 1024)

    assert isinstance(error, LoopLinkError)
    assert error.object_id == "mesh-1"
    assert error.size == 2048
    assert error.limit == 1024
    text = str(error)
    assert "mesh-1" in text
    assert "2048" in text
    assert "1024" in text


def test_viewer_refused_error_keeps_refusals() -> None:
    refusals = [Refusal(object_id="mesh-1", code="too_big", message="Object is too large.")]
    error = ViewerRefusedError(refusals)

    assert isinstance(error, LoopLinkError)
    assert error.refusals == refusals
    assert str(error) != ""


def test_missing_extra_error_keeps_extra_and_gives_install_command() -> None:
    error = MissingExtraError("viewer")

    assert isinstance(error, LoopLinkError)
    assert error.extra == "viewer"
    assert "pip install loop-link[viewer]" in str(error)


def test_refusal_is_hashable_and_frozen() -> None:
    refusal = Refusal(object_id="mesh-1", code="too_big", message="Object is too large.")

    assert hash(refusal) == hash(
        Refusal(object_id="mesh-1", code="too_big", message="Object is too large.")
    )
    assert {refusal}
    with pytest.raises(dataclasses.FrozenInstanceError):
        refusal.code = "other"  # type: ignore[misc]
