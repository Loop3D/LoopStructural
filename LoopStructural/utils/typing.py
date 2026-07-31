import numbers
from typing import TypeVar, Union

T = TypeVar("T")
Array = Union[list[T]]

NumericInput = Union[numbers.Number, Array[numbers.Number]]
