"""
MAPQN Parameters Classes.

Defines parameter structures for various MAPQN analysis algorithms
including linear reduction, MVA-based, and blocking variants.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import List, Optional, Union
from abc import ABC, abstractmethod


class MapqnParameters(ABC):
    """
    Base class for MAPQN model parameters.

    Attributes:
        M: Number of queues.
        N: Population (number of customers).
    """

    @property
    @abstractmethod
    def M(self) -> int:
        """Number of queues."""
        pass

    @property
    @abstractmethod
    def N(self) -> int:
        """Population."""
        pass

    def validate(self) -> None:
        """Validate parameters."""
        if self.M <= 0:
            raise ValueError("M must be positive")
        if self.N <= 0:
            raise ValueError("N must be positive")


@dataclass
class PFParameters(MapqnParameters):
    """
    Parameters for Product Form linear reduction model (no phases).

    Attributes:
        _M: Number of queues.
        _N: Population.
        mu: Service rates array [M].
        r: Routing probability matrix [M x M].
    """
    _M: int
    _N: int
    mu: np.ndarray  # [M]
    r: np.ndarray   # [M x M]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.mu = np.asarray(self.mu, dtype=np.float64).flatten()
        self.r = np.asarray(self.r, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if len(self.mu) != self.M:
            raise ValueError(f"mu array size ({len(self.mu)}) must equal M ({self.M})")
        if self.r.shape != (self.M, self.M):
            raise ValueError(f"r must be {self.M}x{self.M} matrix")
        if np.any(self.mu < 0):
            raise ValueError("Service rates must be non-negative")
        if np.any(self.r < 0):
            raise ValueError("Routing probabilities must be non-negative")

    def q(self, i: int, j: int) -> float:
        """
        Compute q parameter as defined in AMPL model.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).

        Returns:
            Transition rate from i to j.
        """
        return self.r[i, j] * self.mu[i]


@dataclass
class LinearReductionParameters(MapqnParameters):
    """
    Parameters for Linear Reduction models with phases.

    Attributes:
        _M: Number of queues.
        _N: Population.
        K: Number of phases for each queue [M].
        mu: Service rate matrices [M x K[i] x K[i]].
        r: Routing probability matrix [M x M].
        v: Background transition rate matrices [M x K[i] x K[i]].
    """
    _M: int
    _N: int
    K: np.ndarray           # [M] - number of phases per queue
    mu: List[np.ndarray]    # [M] list of K[i] x K[i] matrices
    r: np.ndarray           # [M x M]
    v: List[np.ndarray]     # [M] list of K[i] x K[i] matrices

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.K = np.asarray(self.K, dtype=np.int32).flatten()
        self.r = np.asarray(self.r, dtype=np.float64)
        self.mu = [np.asarray(m, dtype=np.float64) for m in self.mu]
        self.v = [np.asarray(v, dtype=np.float64) for v in self.v]

    def validate(self) -> None:
        super().validate()
        if len(self.K) != self.M:
            raise ValueError(f"K array size must equal M")
        if len(self.mu) != self.M:
            raise ValueError(f"mu array size must equal M")
        if self.r.shape != (self.M, self.M):
            raise ValueError(f"r must be {self.M}x{self.M} matrix")
        if len(self.v) != self.M:
            raise ValueError(f"v array size must equal M")

        for i in range(self.M):
            if self.K[i] <= 0:
                raise ValueError(f"K[{i}] must be positive")
            expected_shape = (self.K[i], self.K[i])
            if self.mu[i].shape != expected_shape:
                raise ValueError(f"mu[{i}] must be {expected_shape} matrix")
            if self.v[i].shape != expected_shape:
                raise ValueError(f"v[{i}] must be {expected_shape} matrix")
            if np.any(self.mu[i] < 0):
                raise ValueError(f"Service rates in mu[{i}] must be non-negative")
            if np.any(self.v[i] < 0):
                raise ValueError(f"Background rates in v[{i}] must be non-negative")

        if np.any(self.r < 0):
            raise ValueError("Routing probabilities must be non-negative")

    def q(self, i: int, j: int, k: int, h: int) -> float:
        """
        Compute q parameter as defined in AMPL model.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).
            k: Source phase index (0-based).
            h: Destination phase index (0-based).

        Returns:
            Transition rate.
        """
        if j != i:
            return self.r[i, j] * self.mu[i][k, h]
        else:
            return self.v[i][k, h] + self.r[i, i] * self.mu[i][k, h]


@dataclass
class MVAVersionParameters(MapqnParameters):
    """
    Parameters for MVA version models (scalar K).

    Attributes:
        _M: Number of queues.
        _N: Population.
        K: Number of levels (scalar).
        muM: Service rates for queues 1..M-1 [M-1].
        muMAP: Service rate matrix for MAP queue [K x K].
        r: Routing probability matrix [M x M].
        v: Level change rate matrix [K x K].
    """
    _M: int
    _N: int
    K: int
    muM: np.ndarray        # [M-1]
    muMAP: np.ndarray      # [K x K]
    r: np.ndarray          # [M x M]
    v: np.ndarray          # [K x K]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.muM = np.asarray(self.muM, dtype=np.float64).flatten()
        self.muMAP = np.asarray(self.muMAP, dtype=np.float64)
        self.r = np.asarray(self.r, dtype=np.float64)
        self.v = np.asarray(self.v, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if self.K <= 0:
            raise ValueError("K must be positive")
        if len(self.muM) != self.M - 1:
            raise ValueError(f"muM array size must equal M-1")
        if self.muMAP.shape != (self.K, self.K):
            raise ValueError(f"muMAP must be {self.K}x{self.K} matrix")
        if self.r.shape != (self.M, self.M):
            raise ValueError(f"r must be {self.M}x{self.M} matrix")
        if self.v.shape != (self.K, self.K):
            raise ValueError(f"v must be {self.K}x{self.K} matrix")
        if np.any(self.muM < 0):
            raise ValueError("Service rates must be non-negative")
        if np.any(self.muMAP < 0):
            raise ValueError("MAP service rates must be non-negative")
        if np.any(self.r < 0):
            raise ValueError("Routing probabilities must be non-negative")
        if np.any(self.v < 0):
            raise ValueError("Level change rates must be non-negative")

    def q(self, i: int, j: int, k: int, h: int) -> float:
        """
        Compute q parameter for MVA model.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).
            k: Source phase/level index (0-based).
            h: Destination phase/level index (0-based).

        Returns:
            Transition rate.
        """
        if i < self.M - 1:
            # Non-MAP queue: scalar service rate
            if k == h:
                return self.r[i, j] * self.muM[i]
            else:
                return 0.0
        else:
            # MAP queue (i == M-1)
            if j < self.M - 1:
                return self.r[self.M - 1, j] * self.muMAP[k, h]
            else:
                if k != h:
                    return self.v[k, h] + self.r[self.M - 1, self.M - 1] * self.muMAP[k, h]
                else:
                    return 0.0


@dataclass
class QRBoundsBasParameters(MapqnParameters):
    """
    Parameters for QR Bounds with Blocking-After-Service (BAS).

    Attributes:
        _M: Number of queues.
        _N: Population.
        MR: Number of independent blocking configurations.
        f: Finite capacity queue index (1-based).
        K: Number of phases for each queue [M].
        F: Capacity for each queue [M].
        MM: Blocking order matrix [MR x M].
        MM1: Blocking order matrix [MR x M].
        ZZ: Non-zeros in independent blocking configurations [MR].
        BB: Blocking state matrix [MR x M].
        mu: Service rate matrices [M] list of K[i] x K[i].
        v: Background transition rate matrices [M] list of K[i] x K[i].
        r: Routing probability matrix [M x M].
    """
    _M: int
    _N: int
    MR: int
    f: int
    K: np.ndarray           # [M]
    F: np.ndarray           # [M]
    MM: np.ndarray          # [MR x M]
    MM1: np.ndarray         # [MR x M]
    ZZ: np.ndarray          # [MR]
    BB: np.ndarray          # [MR x M]
    mu: List[np.ndarray]    # [M] list
    v: List[np.ndarray]     # [M] list
    r: np.ndarray           # [M x M]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    @property
    def ZM(self) -> int:
        """Maximum of ZZ array."""
        return int(np.max(self.ZZ)) if len(self.ZZ) > 0 else 0

    def __post_init__(self):
        self.K = np.asarray(self.K, dtype=np.int32).flatten()
        self.F = np.asarray(self.F, dtype=np.int32).flatten()
        self.MM = np.asarray(self.MM, dtype=np.float64)
        self.MM1 = np.asarray(self.MM1, dtype=np.float64)
        self.ZZ = np.asarray(self.ZZ, dtype=np.int32).flatten()
        self.BB = np.asarray(self.BB, dtype=np.float64)
        self.mu = [np.asarray(m, dtype=np.float64) for m in self.mu]
        self.v = [np.asarray(v, dtype=np.float64) for v in self.v]
        self.r = np.asarray(self.r, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if self.MR <= 0:
            raise ValueError("MR must be positive")
        if not (1 <= self.f <= self.M):
            raise ValueError("f must be between 1 and M")
        if len(self.K) != self.M:
            raise ValueError("K array size must equal M")
        if len(self.F) != self.M:
            raise ValueError("F array size must equal M")

    def q(self, i: int, j: int, k: int, h: int) -> float:
        """Compute q parameter."""
        if j != i:
            return self.r[i, j] * self.mu[i][k, h]
        else:
            return self.v[i][k, h] + self.r[i, i] * self.mu[i][k, h]


@dataclass
class QRBoundsRsrdParameters(MapqnParameters):
    """
    Parameters for QR Bounds with Repetitive-Service Random-Destination (RSRD).

    Attributes:
        _M: Number of queues.
        _N: Population.
        F: Capacity of each queue [M].
        K: Number of phases for each queue [M].
        mu: Service rate matrices [M] list of K[i] x K[i].
        v: Background transition rate matrices [M] list of K[i] x K[i].
        alpha: Load-dependent rates [M x N].
        r: Routing probability matrix [M x M].
    """
    _M: int
    _N: int
    F: np.ndarray           # [M]
    K: np.ndarray           # [M]
    mu: List[np.ndarray]    # [M] list
    v: List[np.ndarray]     # [M] list
    alpha: np.ndarray       # [M x N]
    r: np.ndarray           # [M x M]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.F = np.asarray(self.F, dtype=np.int32).flatten()
        self.K = np.asarray(self.K, dtype=np.int32).flatten()
        self.mu = [np.asarray(m, dtype=np.float64) for m in self.mu]
        self.v = [np.asarray(v, dtype=np.float64) for v in self.v]
        self.alpha = np.asarray(self.alpha, dtype=np.float64)
        self.r = np.asarray(self.r, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if len(self.F) != self.M:
            raise ValueError("F array size must equal M")
        if len(self.K) != self.M:
            raise ValueError("K array size must equal M")
        if len(self.mu) != self.M:
            raise ValueError("mu array size must equal M")
        if len(self.v) != self.M:
            raise ValueError("v array size must equal M")
        if self.alpha.shape[0] != self.M:
            raise ValueError("alpha first dimension must equal M")

    def q(self, i: int, j: int, k: int, h: int, n: int) -> float:
        """
        Compute q parameter with load dependence.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).
            k: Source phase index (0-based).
            h: Destination phase index (0-based).
            n: Population level (1-based in formula, 0 returns 0).

        Returns:
            Transition rate.
        """
        if n == 0:
            return 0.0

        alpha_in = self.alpha[i, n - 1] if n <= self.alpha.shape[1] else 1.0

        if j != i:
            return self.r[i, j] * self.mu[i][k, h] * alpha_in
        else:
            return (self.v[i][k, h] + self.r[i, i] * self.mu[i][k, h]) * alpha_in


@dataclass
class QuadraticDelayParameters(MapqnParameters):
    """
    Parameters for Quadratic Reduction model with delay systems.

    Attributes:
        _M: Number of queues.
        _N: Population.
        K: Number of phases for each queue [M].
        Z: Delay (think time) for delay queue.
        D1: Service demand at queue 1.
        mu: Service rate matrices [M] list of K[i] x K[i].
        v: Background transition rate matrices [M] list of K[i] x K[i].
        alpha: Load-dependent rates [M x N].
        r: Routing probability matrix [M x M].
    """
    _M: int
    _N: int
    K: np.ndarray           # [M]
    Z: float                # Delay (think time)
    D1: float               # Service demand at queue 1
    mu: List[np.ndarray]    # [M] list of K[i] x K[i]
    v: List[np.ndarray]     # [M] list of K[i] x K[i]
    alpha: np.ndarray       # [M x N]
    r: np.ndarray           # [M x M]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.K = np.asarray(self.K, dtype=np.int32).flatten()
        self.mu = [np.asarray(m, dtype=np.float64) for m in self.mu]
        self.v = [np.asarray(v, dtype=np.float64) for v in self.v]
        self.alpha = np.asarray(self.alpha, dtype=np.float64)
        self.r = np.asarray(self.r, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if self.Z <= 0:
            raise ValueError("Z must be positive")
        if self.D1 <= 0:
            raise ValueError("D1 must be positive")
        if len(self.K) != self.M:
            raise ValueError("K array size must equal M")
        if len(self.mu) != self.M:
            raise ValueError("mu array size must equal M")
        if len(self.v) != self.M:
            raise ValueError("v array size must equal M")
        if self.alpha.shape[0] != self.M:
            raise ValueError("alpha first dimension must equal M")
        if self.r.shape != (self.M, self.M):
            raise ValueError(f"r must be {self.M}x{self.M} matrix")

    def q(self, i: int, j: int, k: int, h: int, n: int) -> float:
        """
        Compute q parameter with load dependence.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).
            k: Source phase index (0-based).
            h: Destination phase index (0-based).
            n: Population level (1-based, 0 returns 0).

        Returns:
            Transition rate.
        """
        if n == 0:
            return 0.0

        alpha_in = self.alpha[i, n - 1] if n <= self.alpha.shape[1] else 1.0

        if j != i:
            return self.r[i, j] * self.mu[i][k, h] * alpha_in
        else:
            return (self.v[i][k, h] + self.r[i, i] * self.mu[i][k, h]) * alpha_in


@dataclass
class QuadraticLDParameters(MapqnParameters):
    """
    Parameters for Quadratic Reduction model with load-dependence.

    Attributes:
        _M: Number of queues.
        _N: Population.
        K: Number of phases for each queue [M].
        mu: Service rate matrices [M] list of K[i] x K[i].
        v: Background transition rate matrices [M] list of K[i] x K[i].
        alpha: Load-dependent rates [M x N].
        r: Routing probability matrix [M x M].
    """
    _M: int
    _N: int
    K: np.ndarray           # [M]
    mu: List[np.ndarray]    # [M] list of K[i] x K[i]
    v: List[np.ndarray]     # [M] list of K[i] x K[i]
    alpha: np.ndarray       # [M x N]
    r: np.ndarray           # [M x M]

    @property
    def M(self) -> int:
        return self._M

    @property
    def N(self) -> int:
        return self._N

    def __post_init__(self):
        self.K = np.asarray(self.K, dtype=np.int32).flatten()
        self.mu = [np.asarray(m, dtype=np.float64) for m in self.mu]
        self.v = [np.asarray(v, dtype=np.float64) for v in self.v]
        self.alpha = np.asarray(self.alpha, dtype=np.float64)
        self.r = np.asarray(self.r, dtype=np.float64)

    def validate(self) -> None:
        super().validate()
        if len(self.K) != self.M:
            raise ValueError("K array size must equal M")
        if len(self.mu) != self.M:
            raise ValueError("mu array size must equal M")
        if len(self.v) != self.M:
            raise ValueError("v array size must equal M")
        if self.alpha.shape[0] != self.M:
            raise ValueError("alpha first dimension must equal M")
        if self.r.shape != (self.M, self.M):
            raise ValueError(f"r must be {self.M}x{self.M} matrix")

    def q(self, i: int, j: int, k: int, h: int, n: int) -> float:
        """
        Compute q parameter with load dependence.

        Args:
            i: Source queue index (0-based).
            j: Destination queue index (0-based).
            k: Source phase index (0-based).
            h: Destination phase index (0-based).
            n: Population level (1-based, 0 returns 0).

        Returns:
            Transition rate.
        """
        if n == 0:
            return 0.0

        alpha_in = self.alpha[i, n - 1] if n <= self.alpha.shape[1] else 1.0

        if j != i:
            return self.r[i, j] * self.mu[i][k, h] * alpha_in
        else:
            return (self.v[i][k, h] + self.r[i, i] * self.mu[i][k, h]) * alpha_in


__all__ = [
    'MapqnParameters',
    'PFParameters',
    'LinearReductionParameters',
    'MVAVersionParameters',
    'QRBoundsBasParameters',
    'QRBoundsRsrdParameters',
    'QuadraticDelayParameters',
    'QuadraticLDParameters',
]
