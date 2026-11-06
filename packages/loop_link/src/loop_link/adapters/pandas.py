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
    if values.dtype == object:
        items = list(values)
        if all(isinstance(item, str) for item in items):
            return np.asarray(items, dtype="U")
        raise InvalidInputError(
            f"The column '{name}' has values that are not all strings. "
            "Convert the column to numbers or to strings, then try again."
        )
    if values.dtype.kind in "biuf":
        return values.to_numpy()
    if values.dtype.kind in "US":
        return np.asarray(values.to_numpy(), dtype="U")
    # A string dtype of pandas: all values must be str.
    items = list(values)
    if all(isinstance(item, str) for item in items):
        return np.asarray(items, dtype="U")
    raise InvalidInputError(
        f"The column '{name}' has a type that is not supported. "
        "Convert the column to numbers or to strings, then try again."
    )


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

    if logs is not None:
        logs = logs.reset_index(drop=True)  # Row labels become row positions.
    ids = [str(value) for value in collar[hole_id]]
    known = set(ids)

    def _group(table: pd.DataFrame) -> dict[str, pd.DataFrame]:
        groups: dict[str, pd.DataFrame] = {}
        for key, rows in table.groupby(hole_id, sort=False):
            name = str(key)
            if name not in known:
                raise InvalidInputError(
                    f"The hole id '{name}' is not in the collar table. "
                    "Add the hole to the collar table or remove its rows."
                )
            groups[name] = rows
        return groups

    survey_groups = _group(survey)
    log_groups = _group(logs) if logs is not None else {}

    sign = 1.0 if dip_positive_down else -1.0
    log_columns: list[str] = []
    log_arrays: dict[str, np.ndarray] = {}
    if logs is not None:
        log_columns = [c for c in logs.columns if c not in (hole_id, from_, to)]
        for column in log_columns:
            log_arrays[column] = _column_array(logs[column], str(column))

    holes: dict[str, Drillhole] = {}
    for position, name in enumerate(ids):
        row = collar.iloc[position]
        xyz = [float(row[x]), float(row[y]), float(row[z])]

        rows = survey_groups.get(name)
        if rows is None or len(rows) == 0:
            survey_values = np.zeros((0, 3), dtype=np.float64)
        else:
            rows = rows.sort_values(depth, kind="stable")
            survey_values = np.column_stack(
                [
                    rows[depth].to_numpy(dtype=np.float64),
                    rows[azimuth].to_numpy(dtype=np.float64),
                    sign * rows[dip].to_numpy(dtype=np.float64),
                ]
            )

        interval_values = np.zeros((0, 2), dtype=np.float64)
        interval_properties: dict[str, np.ndarray] = {}
        if logs is not None:
            log_rows = log_groups.get(name)
            if log_rows is not None and len(log_rows) > 0:
                interval_values = np.column_stack(
                    [
                        log_rows[from_].to_numpy(dtype=np.float64),
                        log_rows[to].to_numpy(dtype=np.float64),
                    ]
                )
                index = log_rows.index.to_numpy()
                for column in log_columns:
                    interval_properties[str(column)] = log_arrays[column][index]
            else:
                for column in log_columns:
                    interval_properties[str(column)] = log_arrays[column][:0]

        holes[name] = Drillhole(
            collar=xyz,
            survey=survey_values,
            intervals=interval_values,
            properties={},
            interval_properties=interval_properties,
        )
    return holes
