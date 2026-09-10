"""
Phase-type distribution utilities for FLD solver.

Converts dict-based service distribution representations to matrix-based
phase-type (PH) distributions suitable for fluid ODE analysis.

Phase-type representation:
- D0: Internal transition matrix (sub-generator)
- D1: Completion rate matrix (absorption rates)
- pie: Initial phase probability vector

Supported distributions:
- Exponential: {'rate': mu} -> 1 phase
- Erlang: {'k': k, 'mu': mu} -> k phases
- HyperExponential: {'p': [...], 'mu': [...]} -> n phases
- Cox/PH: {'D0': D0, 'D1': D1} or {'D0': D0, 'pie': pie} -> direct PH

References:
- MATLAB LINE: refreshProcessPhases.m, refreshProcessRepresentations.m
- Phase-type distributions in queueing theory
"""

import numpy as np
from typing import Dict, Tuple, Optional, List, Union, Any
from dataclasses import dataclass


@dataclass
class PhaseTypeRepresentation:
    """Phase-type distribution representation."""
    D0: np.ndarray  # Internal transition matrix (n x n)
    D1: np.ndarray  # Completion rate matrix (n x 1 or n x n)
    pie: np.ndarray  # Initial phase probability vector (n,)
    n_phases: int  # Number of phases

    @property
    def mean(self) -> float:
        """Compute mean of the PH distribution: E[T] = -pie @ inv(D0) @ ones."""
        if self.n_phases == 0:
            return 0.0
        try:
            ones = np.ones(self.n_phases)
            return float(-self.pie @ np.linalg.solve(self.D0, ones))
        except np.linalg.LinAlgError:
            return 0.0

    @property
    def scv(self) -> float:
        """Compute squared coefficient of variation."""
        if self.n_phases == 0:
            return 0.0
        try:
            ones = np.ones(self.n_phases)
            inv_D0 = np.linalg.inv(self.D0)
            m1 = -self.pie @ inv_D0 @ ones
            m2 = 2 * self.pie @ inv_D0 @ inv_D0 @ ones
            var = m2 - m1**2
            return float(var / (m1**2)) if m1 > 0 else 0.0
        except np.linalg.LinAlgError:
            return 1.0


def _map_pie_or_first_phase(D0: np.ndarray, D1: np.ndarray) -> np.ndarray:
    """Arrival-instant phase distribution of (D0, D1), or e_1 if degenerate.

    map_pie is the correct initial vector for any Markovian process: for a PH,
    where D1 = s*alpha has rank one, it returns alpha, so the acyclic case is
    unchanged; for a genuine MAP or MMPP2 it is the embedded distribution seen
    at an arrival, which is not concentrated on phase 1.
    """
    n = int(np.asarray(D0).shape[0])
    try:
        from ....api.mam import map_pie
        pie = np.asarray(map_pie(np.asarray(D0, dtype=float),
                                 np.asarray(D1, dtype=float)), dtype=float).ravel()
        if pie.size == n and np.all(np.isfinite(pie)) and np.all(pie >= -1e-12):
            total = float(np.sum(pie))
            if total > 0:
                return np.maximum(pie, 0.0) / total
    except Exception:
        pass
    pie = np.zeros(n)
    pie[0] = 1.0
    return pie


def convert_dict_to_phase_type(
    proc_dict: Dict[str, Any],
    rate: float = 1.0
) -> PhaseTypeRepresentation:
    """
    Convert dict-based service distribution to phase-type representation.

    Args:
        proc_dict: Dictionary with distribution parameters
        rate: Default rate if not specified

    Returns:
        PhaseTypeRepresentation with D0, D1, pie, n_phases
    """
    if proc_dict is None:
        # Default to exponential with given rate
        return _exponential_ph(rate)

    if not isinstance(proc_dict, dict):
        # Already a list/tuple [D0, D1] - convert to PhaseTypeRepresentation
        if isinstance(proc_dict, (list, tuple)) and len(proc_dict) >= 2:
            D0 = np.asarray(proc_dict[0])
            D1 = np.asarray(proc_dict[1])
            n_phases = D0.shape[0]
            pie = _map_pie_or_first_phase(D0, D1)
            if len(proc_dict) >= 3:
                pie = np.asarray(proc_dict[2]).flatten()
            return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=n_phases)
        # Single matrix - interpret as D0
        if hasattr(proc_dict, 'shape'):
            D0 = np.asarray(proc_dict)
            n_phases = D0.shape[0]
            D1 = -np.sum(D0, axis=1, keepdims=True)
            pie = np.zeros(n_phases)
            pie[0] = 1.0
            return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=n_phases)
        return _exponential_ph(rate)

    # Dict-based format
    if 'D0' in proc_dict:
        # Direct PH representation
        D0 = np.asarray(proc_dict['D0'])
        n_phases = D0.shape[0]
        if 'D1' in proc_dict:
            D1 = np.asarray(proc_dict['D1'])
        else:
            D1 = -np.sum(D0, axis=1, keepdims=True)
        if 'pie' in proc_dict:
            pie = np.asarray(proc_dict['pie']).flatten()
        else:
            pie = np.zeros(n_phases)
            pie[0] = 1.0
        return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=n_phases)

    if 'k' in proc_dict and 'mu' in proc_dict:
        # Erlang(k, mu)
        return _erlang_ph(proc_dict['k'], proc_dict['mu'])

    if 'k' in proc_dict and 'rate' in proc_dict:
        # Erlang with rate parameter (mean = k/rate, so mu = rate)
        return _erlang_ph(proc_dict['k'], proc_dict['rate'])

    if 'rate' in proc_dict:
        # Exponential
        return _exponential_ph(proc_dict['rate'])

    if 'p' in proc_dict and 'mu' in proc_dict:
        # HyperExponential (p/mu format)
        return _hyperexp_ph(proc_dict['p'], proc_dict['mu'])

    if 'probs' in proc_dict and 'rates' in proc_dict:
        # HyperExponential (probs/rates format from HyperExp class)
        return _hyperexp_ph(proc_dict['probs'], proc_dict['rates'])

    if 'lambda' in proc_dict:
        # Alternative exponential notation
        return _exponential_ph(proc_dict['lambda'])

    if 'mean' in proc_dict:
        # Mean-only specification - use exponential
        mean = proc_dict['mean']
        return _exponential_ph(1.0 / mean if mean > 0 else 1.0)

    # Default to exponential
    return _exponential_ph(rate)


def _exponential_ph(rate: float) -> PhaseTypeRepresentation:
    """Create exponential phase-type with given rate."""
    D0 = np.array([[-rate]])
    D1 = np.array([[rate]])
    pie = np.array([1.0])
    return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=1)


def _erlang_ph(k: int, mu: float) -> PhaseTypeRepresentation:
    """
    Create Erlang-k phase-type distribution.

    Erlang(k, mu) has k phases, each with rate mu.
    Mean = k/mu, Variance = k/mu^2, SCV = 1/k

    D0 structure:
        -mu  mu   0   0
         0  -mu  mu   0
         0   0  -mu  mu
         0   0   0  -mu

    D1 structure:
        0
        0
        0
        mu  (only last row has completion)
    """
    k = int(k)
    if k <= 0:
        k = 1
    if mu <= 0:
        mu = 1.0

    D0 = np.zeros((k, k))
    for i in range(k):
        D0[i, i] = -mu
        if i < k - 1:
            D0[i, i + 1] = mu

    D1 = np.zeros((k, 1))
    D1[k - 1, 0] = mu

    pie = np.zeros(k)
    pie[0] = 1.0

    return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=k)


def _hyperexp_ph(p: List[float], mu: List[float]) -> PhaseTypeRepresentation:
    """
    Create HyperExponential phase-type distribution.

    HyperExp with n branches: each branch i has probability p[i] and rate mu[i].
    No transitions between phases - direct absorption.

    D0 structure (diagonal):
        -mu1   0    0
          0  -mu2   0
          0    0  -mu3

    D1 structure:
        mu1
        mu2
        mu3

    pie = [p1, p2, p3]
    """
    p = np.asarray(p).flatten()
    mu = np.asarray(mu).flatten()

    n = len(mu)
    if len(p) != n:
        # Normalize or pad p
        if len(p) < n:
            p = np.concatenate([p, np.zeros(n - len(p))])
        else:
            p = p[:n]

    # Normalize probabilities
    p_sum = np.sum(p)
    if p_sum > 0:
        p = p / p_sum
    else:
        p = np.ones(n) / n

    D0 = np.diag(-mu)
    D1 = mu.reshape(-1, 1)
    pie = p

    return PhaseTypeRepresentation(D0=D0, D1=D1, pie=pie, n_phases=n)


def is_nhpp(sn, i: int, r: int) -> bool:
    """Whether (station i, class r) carries a non-homogeneous Poisson process.

    An NHPP slot of sn.proc holds the rate schedule {breakpoints, rates,
    cyclic}, not a {D0, D1} MAP, so it must never be read as phase-type
    matrices. Detected from sn.procid, as in MATLAB solver_fluid.m.
    """
    from ....constants import ProcessType

    procid = getattr(sn, 'procid', None)
    if procid is None:
        return False
    procid = np.asarray(procid, dtype=object)
    if i >= procid.shape[0] or r >= procid.shape[1]:
        return False
    return procid[i, r] == ProcessType.NHPP


def _procid_is(sn, i: int, r: int, *names) -> bool:
    """Whether sn.procid[i,r] is one of the named ProcessType members."""
    from ....constants import ProcessType

    procid = getattr(sn, 'procid', None)
    if procid is None:
        return False
    procid = np.asarray(procid, dtype=object)
    if i >= procid.shape[0] or r >= procid.shape[1]:
        return False
    return any(procid[i, r] == getattr(ProcessType, n) for n in names)


def is_mapt(sn, i: int, r: int) -> bool:
    """Whether (station i, class r) carries a time-inhomogeneous MAP."""
    return _procid_is(sn, i, r, 'MAPT')


def is_pht(sn, i: int, r: int) -> bool:
    """Whether (station i, class r) carries a time-inhomogeneous phase-type."""
    return _procid_is(sn, i, r, 'PHT')


def is_rate_schedule(sn, i: int, r: int) -> bool:
    """Whether (station i, class r) carries any time-varying process.

    Covers NHPP, MAPt and PHt: their sn.proc slot holds a schedule, not a single
    (D0, D1) pair, so it must never be read as phase-type matrices directly.
    """
    return _procid_is(sn, i, r, 'NHPP', 'MAPT', 'PHT')


def schedule_segments(proc_slot, kind: str):
    """Unpack a MAPt or PHt sn.proc slot into its per-segment matrix pairs.

    The slot is flat, [breakpoints, A_1..A_n, B_1..B_n, cyclic], so that the JAR
    MatrixCell (which cannot nest) carries the same layout; n follows from the
    length. For MAPt the pairs are (D0, D1); for PHt they are (alpha, S) and the
    equivalent MAP pair is (S, s*alpha).

    Returns:
        (breakpoints, pairs, cyclic) with pairs a list of n (A_k, B_k) tuples.
    """
    if proc_slot is None or len(proc_slot) < 4:
        raise ValueError('%s slot of sn.proc must hold at least one segment' % kind)
    n = (len(proc_slot) - 2) // 2
    bp = np.asarray(proc_slot[0], dtype=float).ravel()
    cyclic = bool(proc_slot[-1])
    first = [np.atleast_2d(np.asarray(proc_slot[1 + k], dtype=float)) for k in range(n)]
    second = [np.atleast_2d(np.asarray(proc_slot[1 + n + k], dtype=float)) for k in range(n)]
    return bp, list(zip(first, second)), cyclic


def schedule_map_pair_at(proc_slot, kind: str, seg: int) -> Tuple[np.ndarray, np.ndarray]:
    """The (D0, D1) MAP pair of segment `seg` of a MAPt or PHt slot."""
    _, pairs, _ = schedule_segments(proc_slot, kind)
    a, b = pairs[seg]
    if kind == 'MAPt':
        return a, b
    # PHt: a is a 1-by-h alpha row, b the sub-generator S
    alpha = np.asarray(a, dtype=float).ravel()
    S = np.asarray(b, dtype=float)
    exit_vec = -S.sum(axis=1).reshape(-1, 1)
    return S, exit_vec @ alpha.reshape(1, -1)


def schedule_nominal_pair(proc_slot, kind: str) -> Tuple[np.ndarray, np.ndarray]:
    """Width-weighted time-average (D0, D1) of a MAPt or PHt slot.

    This nominal is the stationary carrier of the phase structure: the fluid ODE
    bakes it into the base rates and expresses each segment as a per-event
    multiplier on it. A convex combination of generators is a generator, so it
    is itself a valid MAP.
    """
    bp, pairs, _ = schedule_segments(proc_slot, kind)
    widths = np.diff(bp)
    total = float(np.sum(widths))
    D0 = None
    D1 = None
    for k in range(len(pairs)):
        d0, d1 = schedule_map_pair_at(proc_slot, kind, k)
        w = float(widths[k]) / total
        D0 = w * d0 if D0 is None else D0 + w * d0
        D1 = w * d1 if D1 is None else D1 + w * d1
    return D0, D1


def prepare_phase_type_structures(sn, segment: Optional[Dict[Tuple[int, int], int]] = None
                                  ) -> Tuple[Dict, Dict, np.ndarray]:
    """
    Prepare phase-type structures from network struct.

    Converts dict-based proc to matrix-based and computes phases.

    Args:
        sn: NetworkStruct with proc, rates

    Returns:
        Tuple of:
        - proc_matrix: Dict[station][class] = [D0, D1]
        - pie_dict: Dict[station][class] = pie vector
        - phases: (M x K) array of phase counts
    """
    M = sn.nstations
    K = sn.nclasses

    phases = np.ones((M, K), dtype=int)
    proc_matrix = {}
    pie_dict = {}

    # Get rates for default exponential
    rates = sn.rates if sn.rates is not None else np.ones((M, K))

    for i in range(M):
        proc_matrix[i] = {}
        pie_dict[i] = {}

        for r in range(K):
            # Get rate for this station-class
            rate = rates[i, r] if i < rates.shape[0] and r < rates.shape[1] else 1.0
            if rate <= 0:
                rate = 1.0

            # Get proc if available
            proc_ir = None
            if hasattr(sn, 'proc') and sn.proc is not None:
                if isinstance(sn.proc, dict):
                    if i in sn.proc and r in sn.proc[i]:
                        proc_ir = sn.proc[i][r]
                elif isinstance(sn.proc, list) and i < len(sn.proc):
                    if sn.proc[i] is not None and r < len(sn.proc[i]):
                        proc_ir = sn.proc[i][r]

            # see _kb/06-solver-catalog.md (Fluid: "Time-varying rates + NHPP
            # transient") for why an NHPP source substitutes a nominal exp process
            if is_nhpp(sn, i, r):
                proc_ir = None

            # MAPt/PHt: substitute a stationary carrier with the SAME phase
            # count, either the segment in force (when the caller asks for one)
            # or the time-averaged nominal. Collapsing to a single exponential,
            # as NHPP does, would discard the phase structure the schedule
            # modulates.
            kind = 'MAPt' if is_mapt(sn, i, r) else ('PHt' if is_pht(sn, i, r) else None)
            if kind is not None:
                seg = None if segment is None else segment.get((i, r))
                if seg is None:
                    D0s, D1s = schedule_nominal_pair(proc_ir, kind)
                else:
                    D0s, D1s = schedule_map_pair_at(proc_ir, kind, int(seg))
                ph = PhaseTypeRepresentation(
                    D0=D0s, D1=D1s, pie=_map_pie_or_first_phase(D0s, D1s),
                    n_phases=int(D0s.shape[0]))
                phases[i, r] = ph.n_phases
                proc_matrix[i][r] = [ph.D0, ph.D1]
                pie_dict[i][r] = ph.pie
                continue

            # Convert to phase-type
            ph = convert_dict_to_phase_type(proc_ir, rate)

            # Store results
            phases[i, r] = ph.n_phases
            proc_matrix[i][r] = [ph.D0, ph.D1]
            pie_dict[i][r] = ph.pie

    return proc_matrix, pie_dict, phases


def extract_mu_phi_from_phase_type(
    proc_matrix: Dict,
    phases: np.ndarray
) -> Tuple[Dict, Dict]:
    """
    Extract Mu and Phi vectors from phase-type matrices.

    Mu{i}{r}(k) = service rate in phase k (diagonal of -D0)
    Phi{i}{r}(k) = completion probability from phase k (D1 row sum / Mu)

    Args:
        proc_matrix: Dict[station][class] = [D0, D1]
        phases: (M x K) phase counts

    Returns:
        Tuple of (Mu, Phi) dicts
    """
    M, K = phases.shape
    Mu = {}
    Phi = {}

    for i in range(M):
        Mu[i] = {}
        Phi[i] = {}

        for r in range(K):
            n_phases = int(phases[i, r])

            if i in proc_matrix and r in proc_matrix[i]:
                D0, D1 = proc_matrix[i][r]
                D0 = np.asarray(D0)
                D1 = np.asarray(D1)

                # Mu = -diag(D0)
                mu_vec = -np.diag(D0)

                # Phi = D1 row sums / Mu (completion probability)
                if D1.ndim == 1:
                    d1_sum = D1
                else:
                    d1_sum = np.sum(D1, axis=1)

                phi_vec = np.zeros(n_phases)
                for k in range(n_phases):
                    if mu_vec[k] > 0:
                        phi_vec[k] = d1_sum[k] / mu_vec[k]
                    else:
                        phi_vec[k] = 1.0

                Mu[i][r] = mu_vec
                Phi[i][r] = phi_vec
            else:
                # Default exponential
                Mu[i][r] = np.array([1.0])
                Phi[i][r] = np.array([1.0])

    return Mu, Phi


def compute_q_indices(phases: np.ndarray) -> np.ndarray:
    """
    Compute starting index in state vector for each station-class.

    q_indices[i, r] = starting index for station i, class r

    Args:
        phases: (M x K) phase counts

    Returns:
        q_indices: (M x K) starting indices (0-based)
    """
    M, K = phases.shape
    q_indices = np.zeros((M, K), dtype=int)

    idx = 0
    for i in range(M):
        for r in range(K):
            q_indices[i, r] = idx
            idx += int(phases[i, r])

    return q_indices


__all__ = [
    'PhaseTypeRepresentation',
    'is_nhpp',
    'is_mapt',
    'is_pht',
    'is_rate_schedule',
    'schedule_segments',
    'schedule_map_pair_at',
    'schedule_nominal_pair',
    'convert_dict_to_phase_type',
    'prepare_phase_type_structures',
    'extract_mu_phi_from_phase_type',
    'compute_q_indices',
]
