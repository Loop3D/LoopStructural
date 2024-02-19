from typing import TypeVar, Union, List
import numbers

T = TypeVar("T")
Array = Union[List[T]]

NumericInput = Union[numbers.Number, Array[numbers.Number]]
