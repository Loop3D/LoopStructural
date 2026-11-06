"""Adapter that builds drillholes from pandas tables."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from loop_link.errors import InvalidInputError, MissingExtraError
from loop_link.objects import Drillhole

if TYPE_CHECKING:
    import pandas as pd


def _column_array(values: pd.Series, name: str) -> np.ndarray:
    """Convert one table column to a numpy array.

    Parameters
    ----------
    values : pandas.Series
        The column values.
    name : str
        The column name. The error text uses it.

    Returns
    -------
    numpy.ndarray
        A ``U`` array for a string column. Other columns keep their dtype.

    Raises
    ------
    InvalidInputError
        If the column is an ``object`` column that holds a value that is not a string.
    """
    if values.dtype.kind != "f" and values.isna().any():
        raise InvalidInputError(
            f"The column '{name}' has missing values. "
            "Fill or remove the missing values, then try again."
        )
    if values.dtype.kind in "biuf" and values.dtype != object:
        return values.to_numpy()
    if values.dtype.kind in "US":
        return np.asarray(values.to_numpy(), dtype="U")
    if values.dtype == object or str(values.dtype) in ("string", "str"):
        if values.dtype == object and values.isna().any():
            raise InvalidInputError(
                f"The column '{name}' has missing values. "
                "Fill or remove the missing values, then try again."
            )
        if pd_infer(values) == "string":
            return np.asarray(values.to_numpy(dtype=object), dtype="U")
        raise InvalidInputError(
            f"The column '{name}' has values that are not all strings. "
            "Convert the column to numbers or to strings, then try again."
        )
    raise InvalidInputError(
        f"The column '{name}' has a type that is not supported. "
        "Convert the column to numbers or to strings, then try again."
    )


def pd_infer(values: pd.Series) -> str:
    """Return the pandas inferred type name of a column."""
    from pandas.api.types import infer_dtype

    return str(infer_dtype(values, skipna=False))


def _require(table: pd.DataFrame, table_name: str, columns: list[str]) -> None:
    """Raise if a table does not have all the named columns."""
    for column in columns:
        if column not in table.columns:
            raise InvalidInputError(
                f"The column '{column}' is not in the {table_name} table. "
                "Check the column names and try again."
            )


def _hole_ids(table: pd.DataFrame, table_name: str, hole_id: str) -> list[str]:
    """Return the hole ids of a table as strings. Raise if an id is missing."""
    column = table[hole_id]
    if column.isna().any():
        raise InvalidInputError(
            f"The {table_name} table has a missing hole id in column '{hole_id}'. "
            "Give each row a hole id."
        )
    return column.astype(str).tolist()


def drillholes_from_pandas(
    collar: pd.DataFrame,
    survey: pd.DataFrame,
    logs: pd.DataFrame | None = None,
    *,
    hole_id: str = "hole_id",
    x: str = "x",
    y: str = "y",
    z: str = "z",
    depth: str = "depth",
    azimuth: str = "azimuth",
    dip: str = "dip",
    from_: str = "from",
    to: str = "to",
    dip_positive_down: bool = True,
) -> dict[str, Drillhole]:
    """Build drillholes from a collar table, a survey table and a log table.

    Parameters
    ----------
    collar : pandas.DataFrame
        One row for each hole, with the hole id and the x, y and z of the collar.
    survey : pandas.DataFrame
        Survey stations, with the hole id, depth, azimuth and dip.
    logs : pandas.DataFrame, optional
        Interval rows, with the hole id, from and to. All other columns become
        interval properties.
    hole_id : str
        Name of the hole id column in all tables.
    x, y, z : str
        Names of the collar coordinate columns.
    depth, azimuth, dip : str
        Names of the survey columns.
    from_, to : str
        Names of the interval columns in ``logs``.
    dip_positive_down : bool
        If True, a dip that points down is positive. If False, the dip sign is
        reversed.

    Returns
    -------
    dict of str to Drillhole
        One drillhole for each row of ``collar``, in the order of ``collar``.
        The keys are the hole ids as strings.

    Raises
    ------
    MissingExtraError
        If pandas is not installed.
    InvalidInputError
        If a hole id in ``survey`` or ``logs`` is not in ``collar``, or if a log
        column has a type that is not supported.
    """
    try:
        import pandas as pd  # noqa: F401
    except ImportError as exc:
        raise MissingExtraError("pandas") from exc

    _require(collar, "collar", [hole_id, x, y, z])
    _require(survey, "survey", [hole_id, depth, azimuth, dip])
    if logs is not None:
        _require(logs, "logs", [hole_id, from_, to])

    ids = _hole_ids(collar, "collar", hole_id)
    if len(set(ids)) != len(ids):
        seen: set[str] = set()
        for name in ids:
            if name in seen:
                raise InvalidInputError(
                    f"The hole id '{name}' is in the collar table more than once. "
                    "Use one collar row for each hole."
                )
            seen.add(name)
    known = set(ids)

    def _group(table: pd.DataFrame, table_name: str, sort_by: str) -> tuple[pd.DataFrame, dict]:
        """Sort a table once and return it with the row positions of each hole."""
        table = table.assign(**{hole_id: _hole_ids(table, table_name, hole_id)})
        table = table.sort_values(sort_by, kind="stable").reset_index(drop=True)
        positions = table.groupby(hole_id, sort=False).indices
        for name in positions:
            if name not in known:
                raise InvalidInputError(
                    f"The hole id '{name}' is not in the collar table. "
                    "Add the hole to the collar table or remove its rows."
                )
        return table, positions

    survey_sorted, survey_pos = _group(survey, "survey", depth)
    sign = 1.0 if dip_positive_down else -1.0
    survey_all = np.column_stack(
        [
            survey_sorted[depth].to_numpy(dtype=np.float64),
            survey_sorted[azimuth].to_numpy(dtype=np.float64),
            sign * survey_sorted[dip].to_numpy(dtype=np.float64),
        ]
    ).reshape(-1, 3)

    log_pos: dict = {}
    intervals_all = np.zeros((0, 2), dtype=np.float64)
    log_arrays: dict[str, np.ndarray] = {}
    if logs is not None:
        logs_sorted, log_pos = _group(logs, "logs", from_)
        intervals_all = np.column_stack(
            [
                logs_sorted[from_].to_numpy(dtype=np.float64),
                logs_sorted[to].to_numpy(dtype=np.float64),
            ]
        ).reshape(-1, 2)
        for column in logs_sorted.columns:
            if column not in (hole_id, from_, to):
                log_arrays[str(column)] = _column_array(logs_sorted[column], str(column))

    xyz = np.column_stack(
        [
            collar[x].to_numpy(dtype=np.float64),
            collar[y].to_numpy(dtype=np.float64),
            collar[z].to_numpy(dtype=np.float64),
        ]
    )
    empty = np.zeros(0, dtype=np.intp)
    holes: dict[str, Drillhole] = {}
    for position, name in enumerate(ids):
        rows = survey_pos.get(name, empty)
        log_rows = log_pos.get(name, empty)
        holes[name] = Drillhole(
            collar=xyz[position],
            survey=survey_all[rows],
            intervals=intervals_all[log_rows],
            properties={},
            interval_properties={col: arr[log_rows] for col, arr in log_arrays.items()},
        )
    return holes
