import numbers
from typing import List, TypeVar, Union

T = TypeVar("T")
Array = Union[List[T]]

NumericInput = Union[numbers.Number, Array[numbers.Number]]
