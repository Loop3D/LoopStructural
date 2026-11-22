from __future__ import annotations

from abc import ABCMeta, abstractmethod

import numpy as np

from ....utils import getLogger

logger = getLogger(__name__)


def smooth_peak(x):
    v = np.zeros(x.shape)
    mask = np.logical_and(x >= -1, x <= 1)
    v[mask] = x[mask] ** 4 - 2 * x[mask] ** 2 + 1
    return v


# smooth_peak is zero for |x| >= 1: a tip at both ends (see FaultDisplacement.tip_ends).
smooth_peak.tip_ends = (1, -1)


class FaultProfileFunction(metaclass=ABCMeta):
    def __init__(self):
        self.lim = [-1, 1]

    @abstractmethod
    def to_dict(self) -> dict:
        pass

    @abstractmethod
    def __call__(self, v: np.ndarray) -> np.ndarray:
        pass

    @property
    def tip_ends(self) -> tuple | None:
        """The ends (+1, -1) of the coordinate beyond which the profile is zero.

        A tip at +1 means the profile is exactly zero for v >= 1. The default
        is None: not declared (``FaultDisplacement.tip_ends`` then tests the
        profile numerically).
        """
        return None

    def plot(self, ax=None):
        if ax is None:
            import matplotlib.pyplot as plt

            _fig, ax = plt.subplots()
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
        """Set the limits of the fault frame coordinate outside of which the
        function value is clamped to the value at the limit

        Parameters
        ----------
        min_x : float
            minimum value of the fault frame coordinate
        max_x : float
            maximum value of the fault frame coordinate
        """
        self.lim = [min_x, max_x]

    def check(self):
        if len(self.B) < 3:
            logger.error("underdetermined")
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

    @property
    def tip_ends(self) -> tuple:
        """The limits at +-1 where the cubic is zero; it is clamped beyond them."""
        self.solve()
        ends = []
        for end, lim in ((1, self.lim[1]), (-1, self.lim[0])):
            value = self._limit_value(lim)
            if lim == end and min(max(value, self.min_v), self.max_v) == 0.0:
                ends.append(end)
        return tuple(ends)

    def _limit_value(self, x: float) -> float:
        """Value at a limit: the value constraint there if there is one, else the cubic.

        The least squares fit leaves about 1e-14 where a constraint asks for 0;
        using the constraint makes a profile that ends at 0 exactly 0 beyond it.
        """
        for cx, cy in self.value_points:
            if cx == x:
                return float(cy)
        return self.w[0] * x**3 + self.w[1] * x**2 + self.w[2] * x + self.w[3]

    def __call__(self, v):
        self.solve()
        eva = self.w[0] * v**3 + self.w[1] * v**2 + self.w[2] * v + self.w[3]
        eva[v >= self.lim[1]] = self._limit_value(self.lim[1])
        eva[v <= self.lim[0]] = self._limit_value(self.lim[0])
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

    tip_ends = ()

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
        data : dict
            Dictionary containing "positive" and "negative" keys, each a dictionary
            of CubicFunction parameters

        Returns
        -------
        Composite
            An initialised composite function given the dictionary parameters
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

    tip_ends = ()

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

    tip_ends = ()

    def to_dict(self) -> dict:
        return {}

    @classmethod
    def from_dict(cls, data: dict) -> Zeros:
        return cls()


class FaultDisplacement:
    def __init__(
        self,
        hw: FaultProfileFunction | None = None,
        fw: FaultProfileFunction | None = None,
        gx: FaultProfileFunction | None = None,
        gy: FaultProfileFunction | None = None,
        gz: FaultProfileFunction | None = None,
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
            logger.info("Gx function none setting to ones")
            self.gx = Ones()
        if self.gy is None:
            logger.info("Gy function none setting to ones")
            self.gy = Ones()
        if self.gz is None:
            logger.info("Gz function none setting to ones")
            self.gz = Ones()

        if self.gx is None:
            raise ValueError("Gx function none can't model fault")
        if self.gy is None:
            raise ValueError("Gy function none can't model fault")
        if self.gz is None:
            raise ValueError("Gz function none can't model fault")

    def __call__(self, gx, gy, gz):

        return self.scale * self.gx(gx) * self.gy(gy) * self.gz(gz)

    @property
    def tip_ends(self) -> dict:
        """The tip of the fault: for gy (1) and gz (2), the ends where the profile is zero.

        A tip at end e of coordinate c means the displacement is exactly zero
        where ``e * g_c >= 1``; ``FaultSegment.evaluate_tip_value`` is built
        from it. A profile declares its ends with an attribute ``tip_ends``
        (the library profiles and ``smooth_peak`` do). For a profile that does
        not, the ends are found by testing that the profile is zero at
        ``|g| = 1, 1.5, 2, 10``.

        Returns
        -------
        dict of int to tuple of int
            ``{1: ends of gy, 2: ends of gz}``; an empty tuple: no tip
        """
        return {1: _profile_tip_ends(self.gy), 2: _profile_tip_ends(self.gz)}

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

    def plot(self, range=(-1, 1), axs: list | None = None):
        try:
            import matplotlib.pyplot as plt

            if axs is None:
                _fig, ax = plt.subplots(1, 3, figsize=(15, 5))
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


def _saturate(distance: np.ndarray) -> np.ndarray:
    """Limit a distance from the fault: the identity up to BLOCK_EXTENSION_EXACT,
    constant beyond BLOCK_EXTENSION_LIMIT, and C2 in between."""
    from ....utils import _fault_block

    exact, limit = _fault_block.BLOCK_EXTENSION_EXACT, _fault_block.BLOCK_EXTENSION_LIMIT
    length = limit - exact
    t = np.clip((distance - exact) / length, 0.0, 1.0)
    # q(0) = 0, q'(0) = 1, q''(0) = 0, q'(1) = 0, q''(1) = 0
    return np.minimum(distance, exact) + length * (t - t**3 + 0.5 * t**4)


def _taylor(function, x: np.ndarray, side: int) -> np.ndarray:
    """Second-order expansion of ``function`` at the fault, from the side ``side``.

    ``d(x0) + d'(x0) (x - x0) + 0.5 d''(x0) (x - x0)**2`` about
    ``x0 = side * BLOCK_SIDE_EPSILON``, the place of a point on the fault in
    the block, with one-sided finite differences of fourth order on six
    points with the step ``BLOCK_EXTENSION_STEP`` (their error is of the
    order of the step to the fourth power times the fifth and sixth
    derivative of the profile; the rounding error of the curvature is about
    1e-9 for a profile of size 1).
    """
    from ....utils import _fault_block

    h = _fault_block.BLOCK_EXTENSION_STEP
    x0 = side * _fault_block.BLOCK_SIDE_EPSILON
    f0, f1, f2, f3, f4, f5 = (function(x0 + side * k * h) for k in range(6))
    slope = side * (-25.0 * f0 + 48.0 * f1 - 36.0 * f2 + 16.0 * f3 - 3.0 * f4) / (12.0 * h)
    curvature = (
        45.0 * f0 - 154.0 * f1 + 214.0 * f2 - 156.0 * f3 + 61.0 * f4 - 10.0 * f5
    ) / (12.0 * h**2)
    return f0 + slope * (x - x0) + 0.5 * curvature * (x - x0) ** 2


def continue_profile(profile, x: np.ndarray, side: int) -> np.ndarray:
    """A gx profile of the side ``side`` continued to the other side of the fault.

    Parameters
    ----------
    profile : callable
        the profile of the block's side: the positive or negative part of a
        :class:`Composite`, or the whole gx profile
    x : np.ndarray
        fault frame coordinate gx of points on the other side
        (``side * x < 0``)
    side : int
        +1 or -1

    Returns
    -------
    np.ndarray
        * :class:`Ones`, :class:`Zeros`: the constant.
        * :class:`CubicFunction`: its own cubic at the signed coordinate
          (without its limits and its floor and ceiling). This is the
          natural continuation, smooth to all orders at the fault.
        * any other callable: the second-order expansion at the fault
          (:func:`_taylor`), so the value, the slope and the curvature are
          continuous (C2) and the third derivative is not.

        Beyond ``BLOCK_EXTENSION_EXACT`` from the fault the coordinate is
        limited (:func:`_saturate`), so the continuation is constant beyond
        ``BLOCK_EXTENSION_LIMIT``.
    """
    from ....utils import _fault_block

    x = np.asarray(x, dtype=float)
    if isinstance(profile, (Ones, Zeros)):
        return np.asarray(profile(x), dtype=float)
    x = -side * _saturate(np.abs(x))
    if isinstance(profile, CubicFunction):
        profile.solve()
        w = profile.w

        def cubic(v):
            return w[0] * v**3 + w[1] * v**2 + w[2] * v + w[3]

        # the cubic is the profile at the fault unless a limit, a floor or a
        # ceiling acts there
        near = side * np.array([_fault_block.BLOCK_SIDE_EPSILON, _fault_block.BLOCK_EXTENSION_STEP])
        if np.allclose(profile(near), cubic(near), rtol=0.0, atol=1e-12):
            return cubic(x)

    def scalar(v):
        return float(np.asarray(profile(np.array([v], dtype=float)), dtype=float)[0])

    return _taylor(scalar, x, side)


def block_displacement(faultfunction, gx, forced_gx, gy, gz, side) -> np.ndarray:
    """The displacement profile of a fault in the block that has the side ``side``.

    On the block's own side (``side * gx >= 0``) this is the profile itself.
    On the other side it is the smooth continuation of the profile of the
    block's side across gx = 0 (``_fault_block.BLOCK_EXTENSION == "smooth"``):

    * a :class:`FaultDisplacement`: its gx profile is continued
      (:func:`continue_profile`); the gy and gz profiles are not changed.
    * any other ``faultfunction(gx, gy, gz)``: the second-order expansion
      in gx at the fault, at the gy and gz of each point.

    The value, the slope and the curvature in gx are continuous at the fault
    (C2), so a field restored with this displacement has a continuous
    gradient across the fault surface. With ``BLOCK_EXTENSION == "mirror"``
    the other side takes the mirror image (the profile at ``side * |gx|``),
    which has a kink where the profile has a slope at the fault.

    Parameters
    ----------
    faultfunction : callable
        ``faultfunction(gx, gy, gz)``
    gx : np.ndarray
        fault frame coordinate of the points
    forced_gx : np.ndarray
        ``force_side(gx, side)``: gx on the block's side
    gy, gz : np.ndarray
    side : int or None
        +1, -1, or None (no block side: the profile at gx)

    Returns
    -------
    np.ndarray
    """
    from ....utils import _fault_block

    out = faultfunction(forced_gx, gy, gz)
    if side not in (-1, 1) or _fault_block.BLOCK_EXTENSION != "smooth":
        return out
    other = side * gx < 0
    if not np.any(other):
        return out
    out = np.array(out, dtype=float, copy=True)
    x, y, z = gx[other], gy[other], gz[other]
    if isinstance(faultfunction, FaultDisplacement):
        profile = faultfunction.gx
        if isinstance(profile, Composite):
            profile = profile.positive if side > 0 else profile.negative
        out[other] = (
            faultfunction.scale
            * continue_profile(profile, x, side)
            * faultfunction.gy(y)
            * faultfunction.gz(z)
        )
    else:
        x = -side * _saturate(np.abs(x))

        def at(v):
            return np.asarray(faultfunction(np.full(x.shape, v), y, z), dtype=float)

        out[other] = _taylor(at, x, side)
    return out


def _profile_tip_ends(profile) -> tuple:
    """The declared tip ends of a profile, or ends found by testing it."""
    ends = getattr(profile, "tip_ends", None)
    if ends is not None:
        return tuple(ends)
    found = []
    for end in (1, -1):
        beyond = end * np.array([1.0, 1.5, 2.0, 10.0])
        if np.all(np.asarray(profile(beyond), dtype=float) == 0.0):
            found.append(end)
    return tuple(found)


class BaseFault:
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


class BaseFault3D:
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
