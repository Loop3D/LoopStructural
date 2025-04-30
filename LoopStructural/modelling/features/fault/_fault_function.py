from __future__ import annotations

from abc import abstractmethod, ABCMeta
from typing import Optional, List
import numpy as np

from ....utils import getLogger

logger = getLogger(__name__)


def smooth_peak(x):
    v = np.zeros(x.shape)
    mask = np.logical_and(x >= -1, x <= 1)
    v[mask] = x[mask] ** 4 - 2 * x[mask] ** 2 + 1
    return v

class FaultProfileFunction(metaclass=ABCMeta):
    def __init__(self):
        self.lim = [-1, 1]
        pass

    @abstractmethod
    def to_dict(self) -> dict:
        pass

    @abstractmethod
    def __call__(self, v: np.ndarray) -> np.ndarray:
        pass

    def plot(self, ax=None):
        if ax is None:
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots()
        x = np.linspace(-1, 1, 100)
        ax.plot(x, self(x), label="ones function")


class CubicFunction(FaultProfileFunction):
    """ """

    def __init__(self):
        """
        Class to represent a cubic function.
        The cubic function is ax**3 + bx**2 + cx + d
        The coefficients a,b,c,d are calculated from the constraints

        """
        super().__init__()
        self.A = []  # np.zeros((4,4))
        self.B = []  # np.zeros((4))
        self.max_v = 999999
        self.min_v = -99999
        self.w = np.zeros(4)
        self.up_to_date = False
        self.value_points = []
        self.gradient_points = []

    def add_cstr(self, x: float, y: float):
        """Add a constraint to the cubic function

        Parameters
        ----------
        x : float
            x value
        y : float
            y value for the function
        """
        self.up_to_date = False
        self.A.append([x**3, x**2, x, 1.0])
        self.B.append(y)
        self.value_points.append([x, y])

    def add_grad(self, x, g):
        """Add a gradient constraint to the cubic function

        Parameters
        ----------
        x : float
            x value
        g : float
            gradient value
        """
        self.up_to_date = False
        self.A.append([3 * x**2, 2 * x, 1.0, 0.0])
        self.B.append(g)
        self.gradient_points.append([x, g])

    def add_max(self, max_v):
        """Adds a ceiling value to the funciton.
        This is used to limit the maximum value returned
        by the function but is not a constraint for the function."""
        self.max_v = max_v

    def add_min(self, min_v):
        """Adds a floor value to the funciton.
        This is used to limit the minimum value returned
        by the function but is not a constraint for the function."""
        self.min_v = min_v

    def set_lim(self, min_x: float, max_x: float):
        """

        Parameters
        ----------
        min_x : _type_
            _description_
        max_x : _type_
            _description_
        """
        self.lim = [min_x, max_x]

    def check(self):
        if len(self.B) < 3:
            print("underdetermined")
            raise ValueError("Underdetermined")

    def solve(self):
        if self.up_to_date:
            return
        self.check()
        A = np.array(self.A)
        B = np.array(self.B)
        ATA = A.T @ A
        ATB = A.T @ B
        self.w = np.linalg.lstsq(ATA, ATB, rcond=None)[0]
        self.up_to_date = True

    def __call__(self, v):
        self.solve()
        eva = self.w[0] * v**3 + self.w[1] * v**2 + self.w[2] * v + self.w[3]
        eva[v > self.lim[1]] = (
            self.w[0] * self.lim[1] ** 3
            + self.w[1] * self.lim[1] ** 2
            + self.w[2] * self.lim[1]
            + self.w[3]
        )
        eva[v < self.lim[0]] = (
            self.w[0] * self.lim[0] ** 3
            + self.w[1] * self.lim[0] ** 2
            + self.w[2] * self.lim[0]
            + self.w[3]
        )
        eva[eva > self.max_v] = self.max_v
        eva[eva < self.min_v] = self.min_v

        return eva

    def to_dict(self) -> dict:
        """Export the function to a dictionary


        Returns
        -------
        dict
            Keys A, B, max_v, min_v, w, up_to_date used to create a new CubicFunction
        """
        return {
            "A": self.A,
            "B": self.B,
            "max_v": self.max_v,
            "min_v": self.min_v,
            "w": self.w.tolist(),
            "value_points": self.value_points,
            "gradient_points": self.gradient_points,
            "up_to_date": self.up_to_date,
        }

    @classmethod
    def from_dict(cls, data: dict) -> CubicFunction:
        """Create a fault profile function from a json dictionary

        Parameters
        ----------
        data : dict
            Dictionary containing A, B, max_v, min_v, w, up_to_date

        Returns
        -------
        CubicFunction
            An initialised function given the dictionary parameters
        """
        instance = cls()
        instance.A = data.get("A", [])
        instance.B = data.get("B", [])
        instance.max_v = data.get("max_v", 999999)
        instance.min_v = data.get("min_v", 999999)
        instance.w = np.array(data.get("w", [0, 0, 0, 0]))
        instance.value_points = data.get("value_points", [])
        instance.gradient_points = data.get("gradient_points", [])
        instance.up_to_date = data.get("up_to_date", False)
        return instance


class Composite(FaultProfileFunction):
    """
    A combination of two profiles for the positive and negative values for a coordinate.
    This is used to model the displacement relative to the fault frame coordinate 0
    """

    def __init__(self, positive: FaultProfileFunction, negative: FaultProfileFunction):
        self.positive = positive
        self.negative = negative

    def __call__(self, v: np.ndarray) -> np.ndarray:
        """calculate the displacement for the input coordinate

        Parameters
        ----------
        v : np.ndarray
            fault frame coordinate between -1 and 1

        Returns
        -------
        np.ndarray
            the displacement for the input coordinate
        """
        v = np.array(v)
        r = np.zeros(v.shape)
        r[v > 0] = self.positive(v[v > 0])
        r[v < 0] = self.negative(v[v < 0])
        return r

    def to_dict(self) -> dict:
        return {
            "positive": self.positive.to_dict(),
            "negative": self.negative.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> Composite:
        """Create a fault profile function from a json dictionary

        Parameters
        ----------
        data : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        positive = CubicFunction.from_dict(data["positive"])
        negative = CubicFunction.from_dict(data["negative"])
        return cls(positive, negative)


class Ones(FaultProfileFunction):
    """
    Returns a fault displacement value of one for the input coordinate
    """

    def __call__(self, v: np.ndarray) -> np.ndarray:
        """calculate the displacement for the input coordinate

        Parameters
        ----------
        v : np.ndarray
            fault frame coordinate between -1 and 1

        Returns
        -------
        np.ndarray
            the displacement for the input coordinate
        """
        v = np.array(v)
        return np.ones(v.shape)

    def to_dict(self) -> dict:
        return {}

    @classmethod
    def from_dict(cls, data: dict) -> Ones:
        return cls()


class Zeros(FaultProfileFunction):
    """
    Returns a fault displacement value of zero for the input coordinate
    """

    def __call__(self, v: np.ndarray) -> np.ndarray:
        """calculate the displacement for the input coordinate

        Parameters
        ----------
        v : np.ndarray
            fault frame coordinate between -1 and 1

        Returns
        -------
        np.ndarray
            the displacement for the input coordinate
        """
        v = np.array(v)
        return np.zeros(v.shape)

    def to_dict(self) -> dict:
        return {}

    @classmethod
    def from_dict(cls, data: dict) -> Zeros:
        return cls()


class FaultDisplacement:
    def __init__(
        self,
        hw: Optional[FaultProfileFunction] = None,
        fw: Optional[FaultProfileFunction] = None,
        gx: Optional[FaultProfileFunction] = None,
        gy: Optional[FaultProfileFunction] = None,
        gz: Optional[FaultProfileFunction] = None,
        scale=0.5,
    ):
        """Function for characterising the displacement of a fault in 3D space
        given the coordinates of the structural frame

        Parameters
        ----------
        hw : Optional[FaultProfileFunction], optional
            hanging wall function, by default None
        fw : Optional[FaultProfileFunction], optional
            footwall function, by default None
        gx : Optional[FaultProfileFunction], optional
            displacement in direction normal to fault surface, by default None
        gy : Optional[FaultProfileFunction], optional
            displacement along fault slip direction, by default None
        gz : Optional[FaultProfileFunction], optional
            direction along fault extent direction, by default None

        """
        self.gx = gx
        if hw is not None and fw is not None:
            self.gx = Composite(hw, fw)
        self.gy = gy
        self.gz = gz
        self.scale = scale
        if self.gx is None:
            print("Gx function none setting to ones")
            self.gx = Ones()
        if self.gy is None:
            print("Gy function none setting to ones")
            self.gy = Ones()
        if self.gz is None:
            print("Gz function none setting to ones")
            self.gz = Ones()

        if self.gx is None:
            raise ValueError("Gx function none can't model fault")
        if self.gy is None:
            raise ValueError("Gy function none can't model fault")
        if self.gz is None:
            raise ValueError("Gz function none can't model fault")

    def __call__(self, gx, gy, gz):

        return self.scale * self.gx(gx) * self.gy(gy) * self.gz(gz)

    def to_dict(self) -> dict:
        return {
            "gx": self.gx.to_dict(),
            "gy": self.gy.to_dict(),
            "gz": self.gz.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict) -> FaultDisplacement:
        gx = CubicFunction.from_dict(data["gx"])
        gy = CubicFunction.from_dict(data["gy"])
        gz = CubicFunction.from_dict(data["gz"])
        return cls(gx=gx, gy=gy, gz=gz)

    def plot(self, range=(-1, 1), axs: Optional[List] = None):
        try:
            import matplotlib.pyplot as plt

            if axs is None:
                fig, ax = plt.subplots(1, 3, figsize=(15, 5))
            for i, (name, f) in enumerate(zip(["gx", "gy", "gz"], [self.gx, self.gy, self.gz])):
                x = np.linspace(range[0], range[1], 100)
                ax[i].plot(x, f(x), label=name)
                ax[i].set_title(name)
                ax[i].set_xlabel("Fault frame coordinate")
                ax[i].set_ylabel("Displacement")
                ax[i].legend()

        except ImportError:
            logger.warning("matplotlib not installed, not plotting")
            return


class BaseFault(object):
    """ """

    hw = CubicFunction()
    hw.add_cstr(0, 1)
    hw.add_grad(0, 0)
    hw.add_cstr(1, 0)
    # hw.add_cstr(1,1)

    hw.add_grad(1, 0)
    hw.set_lim(0, 1)
    fw = CubicFunction()
    fw.add_cstr(0, -1)
    fw.add_grad(0, 0)
    fw.add_cstr(-1, 0)
    fw.add_grad(-1, 0)
    fw.set_lim(-1, 0)
    # gyf = CubicFunction()
    # gyf.add_cstr(-1, 0)
    # gyf.add_cstr(1, 0)
    # gyf.add_cstr(-0.2, 1)
    # gyf.add_cstr(0.2, 1)
    # gyf.add_grad(0, 0)
    # gyf.add_min(-1)
    # gyf.add_max(1)
    gyf = Ones()
    gzf = smooth_peak
    gxf = Composite(hw, fw)
    fault_displacement = FaultDisplacement(gx=gxf, gy=gyf, gz=gzf)


class BaseFault3D(object):
    """ """

    hw = CubicFunction()
    hw.add_cstr(0, 1)
    hw.add_grad(0, 0)
    hw.add_cstr(1, 0)
    # hw.add_cstr(1,1)

    hw.add_grad(1, 0)
    hw.add_max(1)
    fw = CubicFunction()
    fw.add_cstr(0, -1)
    fw.add_grad(0, 0)
    fw.add_cstr(-1, 0)
    fw.add_grad(-1, 0)
    fw.add_min(-1)

    gyf = smooth_peak
    # CubicFunction()
    # gyf.add_cstr(-1, 0)
    # gyf.add_cstr(1, 0)
    # gyf.add_cstr(-0.2, 1)
    # gyf.add_cstr(0.2, 1)
    # gyf.add_grad(0, 0)
    # gyf.add_min(-1)
    # gyf.add_max(1)
    # gyf = Ones()
    gzf = smooth_peak
    gxf = Composite(hw, fw)
    fault_displacement = FaultDisplacement(gx=gxf, gy=gyf, gz=gzf)
