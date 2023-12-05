from typing import TypeVar, Union, List
import numbers
import numpy

T = TypeVar("T")
Array = Union[List[T], numpy.ndarray[T]]

NumericInput = Union[numbers.Number, Array[numbers.Number]]
