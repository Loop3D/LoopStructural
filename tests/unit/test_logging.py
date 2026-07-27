"""Tests for the structured logging/timing infrastructure (ROADMAP.md Stage 1b)."""

import logging
import time

import pytest

import LoopStructural
from LoopStructural.utils import (
    FileSink,
    LogSink,
    SqliteSink,
    StreamSink,
    add_sink,
    getLogger,
    get_levels,
    remove_sink,
    timed,
    timed_stage,
)


@pytest.fixture(autouse=True)
def _clean_extra_sinks():
    """`add_sink` mutates package-level state; keep tests isolated from each other."""
    before = list(LoopStructural._extra_sinks)
    yield
    for handler in list(LoopStructural._extra_sinks):
        if handler not in before:
            remove_sink(handler)


def test_get_levels_contains_expected_keys():
    levels = get_levels()
    assert levels["info"] == logging.INFO
    assert levels["warning"] == logging.WARNING
    assert levels["error"] == logging.ERROR
    assert levels["debug"] == logging.DEBUG


def test_getlogger_returns_stdlib_logger_and_registers_it():
    logger = getLogger("loopstructural.test.plain")
    assert isinstance(logger, logging.Logger)
    assert logger.propagate is False
    assert LoopStructural.loggers["loopstructural.test.plain"] is logger


def test_add_sink_with_plain_callable_receives_records():
    received = []
    handler = add_sink(lambda record: received.append(record.getMessage()))
    logger = getLogger("loopstructural.test.callable_sink")
    logger.setLevel(logging.INFO)
    logger.warning("hello from callable sink")

    assert "hello from callable sink" in received
    remove_sink(handler)


def test_remove_sink_detaches_handler():
    received = []
    handler = add_sink(lambda record: received.append(record.getMessage()))
    logger = getLogger("loopstructural.test.remove_sink")
    remove_sink(handler)
    logger.warning("should not be captured")

    assert received == []


def test_add_sink_attaches_to_loggers_created_afterwards():
    received = []
    handler = add_sink(lambda record: received.append(record.getMessage()))
    # created *after* add_sink -- should still pick up the sink automatically.
    logger = getLogger("loopstructural.test.late_logger")
    logger.warning("late logger message")

    assert "late logger message" in received
    remove_sink(handler)


class _ListSink(LogSink):
    """Minimal LogSink subclass used to test the ABC extension point."""

    def __init__(self):
        self.records = []

    def emit(self, record):
        self.records.append(record)


def test_logsink_subclass_extension_point():
    sink = _ListSink()
    handler = add_sink(sink)
    logger = getLogger("loopstructural.test.logsink_subclass")
    logger.warning("via subclass")

    assert len(sink.records) == 1
    assert sink.records[0].getMessage() == "via subclass"
    remove_sink(handler)


def test_stream_sink_writes_to_given_stream():
    import io

    stream = io.StringIO()
    sink = StreamSink(stream, level=logging.INFO)
    handler = add_sink(sink)
    logger = getLogger("loopstructural.test.stream_sink")
    logger.setLevel(logging.INFO)
    logger.info("streamed message")

    assert "streamed message" in stream.getvalue()
    remove_sink(handler)


def test_file_sink_writes_to_file(tmp_path):
    path = tmp_path / "loop.log"
    sink = FileSink(path, level=logging.INFO)
    handler = add_sink(sink)
    logger = getLogger("loopstructural.test.file_sink")
    logger.setLevel(logging.INFO)
    logger.info("logged to file")
    handler.flush()

    assert "logged to file" in path.read_text()
    remove_sink(handler)


def test_file_sink_creates_parent_directories(tmp_path):
    path = tmp_path / "nested" / "dir" / "loop.log"
    FileSink(path)
    assert path.parent.is_dir()


def test_sqlite_sink_records_are_queryable(tmp_path):
    sink = SqliteSink(tmp_path / "loop.sqlite")
    handler = add_sink(sink)
    logger = getLogger("loopstructural.test.sqlite_sink")
    logger.setLevel(logging.INFO)

    with timed_stage(logger, "example_stage"):
        time.sleep(0.001)

    rows = sink.query(stage="example_stage")
    assert len(rows) == 2  # start + end
    events = {row["event"] for row in rows}
    assert events == {"start", "end"}

    end_row = next(row for row in rows if row["event"] == "end")
    assert end_row["duration_s"] > 0
    assert end_row["logger_name"] == "loopstructural.test.sqlite_sink"

    remove_sink(handler)


def test_sqlite_sink_query_filters_by_run_id(tmp_path):
    sink = SqliteSink(tmp_path / "loop.sqlite")
    handler = add_sink(sink)
    logger = getLogger("loopstructural.test.sqlite_run_id")
    logger.setLevel(logging.INFO)

    with timed_stage(logger, "stage_a", run_id="run-1"):
        pass
    with timed_stage(logger, "stage_b", run_id="run-2"):
        pass

    assert len(sink.query(run_id="run-1")) == 2
    assert len(sink.query(run_id="run-2")) == 2
    assert len(sink.query(run_id="run-1", stage="stage_b")) == 0

    remove_sink(handler)


def test_timed_stage_logs_start_and_end_even_on_exception():
    received = []
    handler = add_sink(lambda record: received.append(getattr(record, "event", None)))
    logger = getLogger("loopstructural.test.timed_stage_exception")
    logger.setLevel(logging.INFO)

    with pytest.raises(ValueError):
        with timed_stage(logger, "failing_stage"):
            raise ValueError("boom")

    assert received == ["start", "end"]
    remove_sink(handler)


def test_timed_decorator_times_a_function_call(tmp_path):
    sink = SqliteSink(tmp_path / "loop.sqlite")
    handler = add_sink(sink)

    @timed("decorated_stage", logger=getLogger("loopstructural.test.timed_decorator"))
    def work(x):
        return x * 2

    getLogger("loopstructural.test.timed_decorator").setLevel(logging.INFO)
    assert work(21) == 42

    rows = sink.query(stage="decorated_stage")
    assert len(rows) == 2
    remove_sink(handler)


def test_log_to_console_and_log_to_file_still_work(tmp_path):
    # Backward-compat smoke test for the pre-existing public functions.
    # log_to_file/log_to_console only reconfigure already-registered loggers
    # (pre-existing behavior, unchanged by Stage 1b) so register first.
    from LoopStructural.utils import log_to_console, log_to_file

    logger = getLogger("loopstructural.test.compat_file")
    logfile = tmp_path / "compat.log"
    log_to_file(str(logfile), level="info")
    logger.info("compat message")
    assert "compat message" in logfile.read_text()

    log_to_console(level="warning")
