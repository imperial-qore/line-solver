"""
Symbolic export of the mean-field ODE system integrated by SolverFLD.

Mirrors the MATLAB implementation in solver_fluid_symodes.m and
SolverFLD/exportODEs.m so that all codebases emit the same LaTeX document
for a given model. Two representations are produced:

form = 'W': dx/dt = W'*theta(x) + lambda   (methods: default, matrix, pnorm)
form = 'J': dx/dt = J*r(x)                 (methods: closing, statedep, softmin)

No JPype or JVM dependencies.
"""

import numpy as np

from ...api.sn import SchedStrategy
from ...api.mc.dtmc import dtmc_stochcomp
from ...api.mam.map_analysis import map_pie

# LaTeX text ids aligned with MATLAB SchedStrategy.toText
_SCHED_TEXT = {
    'INF': 'inf', 'FCFS': 'fcfs', 'FCFSPR': 'fcfspr', 'FCFSPI': 'fcfspi',
    'LCFS': 'lcfs', 'LCFSPR': 'lcfspr', 'LCFSPI': 'lcfspi', 'POLLING': 'polling',
    'SIRO': 'siro', 'SJF': 'sjf', 'LJF': 'ljf', 'PS': 'ps', 'DPS': 'dps',
    'GPS': 'gps', 'SEPT': 'sept', 'LEPT': 'lept', 'SRPT': 'srpt',
    'SRPTPRIO': 'srptprio', 'HOL': 'hol', 'FORK': 'fork', 'EXT': 'ext',
    'REF': 'ref', 'PSPRIO': 'psprio', 'DPSPRIO': 'dpsprio', 'GPSPRIO': 'gpsprio',
    'LCFSPRIO': 'lcfsprio', 'LCFSPRPRIO': 'lcfsprprio', 'LCFSPIPRIO': 'lcfspiprio',
}


class SymODEs:
    """Structural description of the exported mean-field ODE system.

    State variable, station and class indices are 0-based; they are
    rendered 1-based in the LaTeX output.
    """

    def __init__(self):
        self.form = None
        self.method = None
        self.nstates = 0
        self.stateStation = None
        self.stateClass = None
        self.statePhase = None
        self.stationNames = []
        self.classNames = []
        self.schedNames = []
        self.sched = None
        self.S = None
        # form W
        self.W = None
        self.Alambda = None
        self.isSource = None
        self.smoothing = None
        self.pstar = None
        self.keep = None
        # form J
        self.nevents = 0
        self.eventFrom = None
        self.eventTo = None
        self.coeff = None
        self.eventVar = None
        self.factorType = None
        self.factorStation = None
        self.factorClass = None
        self.factorOthers = None
        self.dpsw = None
        self.fcfsPhaseW = None
        self.alpha = None
        # initial condition
        self.x0 = None


class _PhStructs:
    """Phase-type view of the network struct, as prepared for the FLD
    engine: proc[i][r] = [D0, D1], pie[i][r], mu[i][r], phi[i][r], phases."""

    def __init__(self, sn):
        from .utils.phase_type import (prepare_phase_type_structures,
                                       extract_mu_phi_from_phase_type)
        self.proc, self.pie, self.phases = prepare_phase_type_structures(sn)
        self.mu, self.phi = extract_mu_phi_from_phase_type(self.proc, self.phases)

    def enabled(self, sn, i, r):
        return not np.isnan(sn.rates[i, r])

    def proc_full(self, i, r):
        """(D0, D1) MAP pair with D1 as a full matrix. The prepared proc
        stores the exit-rate column; the full D1 closes the process with
        the entry vector, D1 = exit_rates * pie, as in MATLAB sn.proc."""
        D0 = np.asarray(self.proc[i][r][0], dtype=float)
        D1 = np.asarray(self.proc[i][r][1], dtype=float)
        if D1.shape != D0.shape:
            pie = np.asarray(self.pie[i][r], dtype=float).flatten()
            D1 = np.outer(D1.flatten(), pie)
        return D0, D1


def solver_fluid_symodes(sn, options):
    """Build a symbolic description of the mean-field ODE system.

    Args:
        sn: NetworkStruct
        options: SolverFLDOptions (method, pstar)

    Returns:
        SymODEs instance
    """
    M = sn.nstations
    K = sn.nclasses

    method = str(options.method).replace('fluid.', '')
    if method in ('default', 'matrix', 'pnorm'):
        form = 'W'
        if method == 'default':
            method = 'matrix'
    elif method in ('closing', 'statedep', 'softmin'):
        form = 'J'
    else:
        raise ValueError(
            "Symbolic ODE export is unsupported for method '%s'. "
            "Supported methods: default, matrix, pnorm, closing, statedep, softmin." % method)

    sys = SymODEs()
    sys.form = form
    sys.method = method
    sys.sched = [sn.sched.get(i) for i in range(M)]
    station_to_node = np.asarray(sn.stationToNode, dtype=float).flatten()
    for i in range(M):
        sys.stationNames.append(sn.nodenames[int(station_to_node[i])])
        sched_name = sys.sched[i].name if hasattr(sys.sched[i], 'name') else str(sys.sched[i])
        sys.schedNames.append(_SCHED_TEXT.get(sched_name, sched_name.lower()))
    for r in range(K):
        sys.classNames.append(sn.classnames[r])

    ph = _PhStructs(sn)
    if form == 'W':
        _build_wform(sys, sn, ph, options, M, K)
    else:
        _build_jform(sys, sn, ph, M, K, method)
    _build_x0(sys, sn, M, K)
    return sys


def _build_wform(sys, sn, ph, options, M, K):
    """Mirror the ODE construction of the matrix method (Ruuskanen et al.,
    PEVA 151 (2021)): dx/dt = W'*theta(x) + Alambda."""
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    S = np.zeros(M)
    for i in range(M):
        si = float(np.asarray(sn.nservers).flatten()[i])
        S[i] = np.sum(njobs) if np.isinf(si) else si

    phases = np.asarray(ph.phases, dtype=int)

    # station-to-station routing matrix via stochastic complementation
    station_to_stateful = np.asarray(sn.stationToStateful, dtype=float).flatten()
    station_indices = []
    for ist in range(M):
        isf = int(station_to_stateful[ist])
        for r in range(K):
            station_indices.append(isf * K + r)
    P = dtmc_stochcomp(np.asarray(sn.rt, dtype=float), np.asarray(station_indices))

    # remove Sink->Source feedback routing for open classes
    for src_ist in range(M):
        if sys.sched[src_ist] == SchedStrategy.EXT:
            for r in range(K):
                if not np.isnan(sn.rates[src_ist, r]) and sn.rates[src_ist, r] > 0:
                    src_col = src_ist * K + r
                    for from_ist in range(M):
                        if from_ist != src_ist:
                            for from_r in range(K):
                                P[from_ist * K + from_r, src_col] = 0.0

    # W = Psi + B*P*A'
    blocks_psi = []
    blocks_a = []
    blocks_b = []
    for ist in range(M):
        for r in range(K):
            if phases[ist, r] == 0:
                blocks_psi.append(np.zeros((1, 1)))
                blocks_a.append(np.full((1, 1), np.nan))
                blocks_b.append(np.zeros((1, 1)))
            else:
                D0, D1 = ph.proc_full(ist, r)
                pie = np.asarray(ph.pie[ist][r], dtype=float).flatten()
                blocks_psi.append(D0)
                blocks_a.append(pie.reshape(-1, 1))
                blocks_b.append(np.sum(D1, axis=1).reshape(-1, 1))
    psi = _blkdiag(blocks_psi)
    A = _blkdiag(blocks_a)
    B = _blkdiag(blocks_b)
    W = psi + B @ P @ A.T

    # exogenous arrival rates into queue phases
    source_arrivals = np.zeros((M, K))
    for src_ist in range(M):
        if sys.sched[src_ist] == SchedStrategy.EXT:
            for r in range(K):
                if not np.isnan(sn.rates[src_ist, r]) and sn.rates[src_ist, r] > 0:
                    source_arrivals[src_ist, r] = sn.rates[src_ist, r]
    total_states = W.shape[0]
    alambda_full = np.zeros(total_states)
    state = 0
    for ist in range(M):
        for r in range(K):
            nphases = int(phases[ist, r])
            if nphases > 0:
                if sys.sched[ist] == SchedStrategy.EXT:
                    state += nphases
                else:
                    arrival_rate_to_queue = 0.0
                    for src_ist in range(M):
                        if source_arrivals[src_ist, r] > 0:
                            arrival_rate_to_queue += source_arrivals[src_ist, r] * P[src_ist * K + r, ist * K + r]
                    if arrival_rate_to_queue > 0:
                        pie = np.asarray(ph.pie[ist][r], dtype=float).flatten()
                        for k in range(nphases):
                            alambda_full[state] = pie[k] * arrival_rate_to_queue
                            state += 1
                    else:
                        state += nphases
            else:
                state += 1

    # state metadata over the same enumeration, then keep-filter
    st_station = np.zeros(total_states, dtype=int)
    st_class = np.zeros(total_states, dtype=int)
    st_phase = np.zeros(total_states, dtype=int)
    state = 0
    for ist in range(M):
        for r in range(K):
            nphases = int(phases[ist, r])
            if nphases == 0:
                st_station[state] = ist
                st_class[state] = r
                st_phase[state] = 0  # placeholder, removed by keep filter
                state += 1
            else:
                for k in range(nphases):
                    st_station[state] = ist
                    st_class[state] = r
                    st_phase[state] = k + 1
                    state += 1

    keep = np.where(~np.isnan(np.sum(W, axis=0)))[0]
    sys.keep = keep
    sys.W = W[np.ix_(keep, keep)]
    sys.Alambda = alambda_full[keep]
    sys.stateStation = st_station[keep]
    sys.stateClass = st_class[keep]
    sys.statePhase = st_phase[keep]
    sys.nstates = len(keep)
    sys.S = S
    sys.isSource = np.array([sys.sched[s] == SchedStrategy.EXT for s in sys.stateStation])

    # smoothing selection mirrors the pstar gate of the matrix method
    pstar = getattr(options, 'pstar', None)
    if pstar is not None:
        pstar_arr = np.asarray(pstar, dtype=float).flatten()
        if pstar_arr.size == 1:
            pstar_arr = np.full(M, pstar_arr[0])
        sys.smoothing = 'pnorm'
        sys.pstar = pstar_arr
    else:
        sys.smoothing = 'min'
        sys.pstar = None


def _blkdiag(blocks):
    rows = sum(b.shape[0] for b in blocks)
    cols = sum(b.shape[1] for b in blocks)
    out = np.zeros((rows, cols))
    r0 = 0
    c0 = 0
    for b in blocks:
        out[r0:r0 + b.shape[0], c0:c0 + b.shape[1]] = b
        r0 += b.shape[0]
        c0 += b.shape[1]
    return out


def _build_jform(sys, sn, ph, M, K, method):
    """Mirror the event enumeration of the closing/statedep/softmin ODE
    methods: dx/dt = J*r(x)."""
    S = np.zeros(M)
    for i in range(M):
        si = float(np.asarray(sn.nservers).flatten()[i])
        S[i] = sn.nclosedjobs if np.isinf(si) else si

    weighted = method in ('statedep', 'softmin')
    if weighted and any(s == SchedStrategy.EXT for s in sys.sched):
        raise ValueError(
            "The '%s' ODE method does not support open models, so their ODE "
            "system cannot be exported. Use the 'matrix' or 'closing' method instead." % method)

    # state indexing as in solver_fluid_odes
    Mu = [[None] * K for _ in range(M)]
    Phi = [[None] * K for _ in range(M)]
    q_indices = np.zeros((M, K), dtype=int)
    Kic = np.zeros((M, K), dtype=int)
    enabled = np.zeros((M, K), dtype=bool)
    cs = 0
    for i in range(M):
        for c in range(K):
            q_indices[i, c] = cs
            if not ph.enabled(sn, i, c):
                Kic[i, c] = 0
            else:
                Mu[i][c] = np.asarray(ph.mu[i][c], dtype=float).flatten()
                Phi[i][c] = np.asarray(ph.phi[i][c], dtype=float).flatten()
                Kic[i, c] = len(Mu[i][c])
                enabled[i, c] = True
            cs += Kic[i, c]
    nstates = cs

    sys.stateStation = np.zeros(nstates, dtype=int)
    sys.stateClass = np.zeros(nstates, dtype=int)
    sys.statePhase = np.zeros(nstates, dtype=int)
    for i in range(M):
        for c in range(K):
            for k in range(Kic[i, c]):
                sys.stateStation[q_indices[i, c] + k] = i
                sys.stateClass[q_indices[i, c] + k] = c
                sys.statePhase[q_indices[i, c] + k] = k + 1

    # normalized DPS weights
    dpsw = np.ones((M, K))
    for i in range(M):
        if sys.sched[i] == SchedStrategy.DPS:
            dpsw[i, :] = np.asarray(sn.schedparam, dtype=float)[i, :]
            dpsw[i, :] = dpsw[i, :] / np.sum(dpsw[i, :])

    # FCFS phase weights used by statedep/softmin: w = -1/D0(k,k)
    fcfs_phase_w = np.zeros(nstates)
    if weighted:
        for i in range(M):
            if sys.sched[i] == SchedStrategy.FCFS:
                for c in range(K):
                    if enabled[i, c]:
                        D0 = np.asarray(ph.proc[i][c][0], dtype=float)
                        for k in range(Kic[i, c]):
                            fcfs_phase_w[q_indices[i, c] + k] = -1.0 / D0[k, k]

    handled = (SchedStrategy.INF, SchedStrategy.EXT, SchedStrategy.PS,
               SchedStrategy.FCFS, SchedStrategy.DPS)

    # see _kb/06-solver-catalog.md (Fluid: "sn.rt omits pseudo-closed Sink -> Source feedback")
    rt = np.asarray(sn.rt, dtype=float).copy()
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    for c in range(K):
        if np.isinf(njobs[c]):
            src = -1
            for i in range(M):
                if sys.sched[i] == SchedStrategy.EXT and not np.isnan(sn.rates[i, c]) and sn.rates[i, c] > 0:
                    src = i
                    break
            if src >= 0:
                for i in range(M):
                    if i == src or not enabled[i, c]:
                        continue
                    row = i * K + c
                    deficit = 1.0 - np.sum(rt[row, :])
                    if deficit > 1e-12:
                        rt[row, src * K + c] += deficit

    ev_from = []
    ev_to = []
    ev_coeff = []
    ev_var = []
    ev_type = []
    ev_station = []
    ev_class = []
    ev_others = []

    def add_event(from_idx, to_idx, base, schedi, i, c, ki):
        ftype = None
        others = None
        coeff = base
        if method == 'closing':
            if schedi == SchedStrategy.INF:
                ftype = 'lin'
            elif schedi == SchedStrategy.EXT:
                if ki == 0:
                    ftype = 'ext1'
                    others = [q_indices[i, c] + u for u in range(1, Kic[i, c])]
                else:
                    ftype = 'lin'
            elif schedi in (SchedStrategy.PS, SchedStrategy.FCFS):
                ftype = 'min'
            elif schedi == SchedStrategy.DPS:
                # the share w_ir*x/ntilde_i of the capacity min(n_i,S_i), as in
                # the closing rate factors: no additive seed, not the full S_i
                ftype = 'dpsmin'
                coeff = coeff * dpsw[i, c]
            else:
                # strategies without a case in the closing rates keep rates = x
                ftype = 'lin'
        else:  # statedep, softmin
            if schedi == SchedStrategy.INF:
                ftype = 'lin'
            elif schedi == SchedStrategy.PS:
                ftype = 'min'
            elif schedi == SchedStrategy.FCFS:
                ftype = 'fcfsws' if method == 'softmin' else 'fcfsw'
                coeff = coeff * fcfs_phase_w[q_indices[i, c] + ki]
            else:  # DPS
                ftype = 'dpspw'
        ev_from.append(from_idx)
        ev_to.append(to_idx)
        ev_coeff.append(coeff)
        ev_var.append(q_indices[i, c] + ki)
        ev_type.append(ftype)
        ev_station.append(i)
        ev_class.append(c)
        ev_others.append(others)

    for i in range(M):
        for c in range(K):
            if enabled[i, c]:
                xic = q_indices[i, c]
                for j in range(M):
                    for l in range(K):
                        if rt[i * K + c, j * K + l] > 0:
                            if not ph.enabled(sn, j, l):
                                pie = np.array([1.0])
                            else:
                                D0jl, D1jl = ph.proc_full(j, l)
                                pie = np.asarray(map_pie(D0jl, D1jl)).flatten()
                            xjl = q_indices[j, l]
                            for ki in range(Kic[i, c]):
                                for kj in range(Kic[j, l]):
                                    if weighted:
                                        if sys.sched[i] == SchedStrategy.INF and j == i:
                                            continue  # self-loop departures skipped at INF stations
                                        if sys.sched[i] not in handled:
                                            continue
                                    base = Phi[i][c][ki] * Mu[i][c][ki] * rt[i * K + c, j * K + l] * pie[kj]
                                    if base > 0:
                                        add_event(xic + ki, xjl + kj, base, sys.sched[i], i, c, ki)
    for i in range(M):
        for c in range(K):
            if enabled[i, c]:
                if weighted and sys.sched[i] not in handled:
                    continue
                xic = q_indices[i, c]
                D0 = np.asarray(ph.proc[i][c][0], dtype=float)
                for ki in range(Kic[i, c] - 1):
                    for kip in range(Kic[i, c]):
                        if ki != kip:
                            base = D0[ki, kip]
                            if base > 0:
                                add_event(xic + ki, xic + kip, base, sys.sched[i], i, c, ki)

    sys.nevents = len(ev_coeff)
    sys.eventFrom = np.array(ev_from, dtype=int)
    sys.eventTo = np.array(ev_to, dtype=int)
    sys.coeff = np.array(ev_coeff, dtype=float)
    sys.eventVar = np.array(ev_var, dtype=int)
    sys.factorType = ev_type
    sys.factorStation = np.array(ev_station, dtype=int)
    sys.factorClass = np.array(ev_class, dtype=int)
    sys.factorOthers = ev_others
    sys.nstates = nstates
    sys.S = S
    sys.dpsw = dpsw
    sys.fcfsPhaseW = fcfs_phase_w
    if method == 'softmin':
        sys.alpha = 20.0  # softmin parameter, as in the softmin ODE method


def _class_mass_at_station(sn, ist, r, K):
    """Jobs of class r initially at station ist, from the model state
    marginals when available, else from the reference station default."""
    state = getattr(sn, 'state', None)
    if state is not None and len(state) > 0:
        station_to_stateful = np.asarray(sn.stationToStateful, dtype=float).flatten()
        isf = int(station_to_stateful[ist])
        state_vec = None
        if isinstance(state, dict):
            if isf in state and state[isf] is not None:
                state_vec = np.asarray(state[isf], dtype=float).flatten()
        elif isf < len(state) and state[isf] is not None:
            state_vec = np.asarray(state[isf], dtype=float).flatten()
        if state_vec is not None and r < len(state_vec):
            val = float(state_vec[r])
            return 0.0 if np.isnan(val) or np.isinf(val) else val
    refstat = np.asarray(sn.refstat, dtype=float).flatten()
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    if int(refstat[r]) == ist and np.isfinite(njobs[r]):
        return float(njobs[r])
    return 0.0


def _build_x0(sys, sn, M, K):
    """Initial condition in the exported state space: phase-1 placement of
    the initial per-class station masses, as in solver_fluid_initsol."""
    x0 = np.zeros(sys.nstates)
    for s in range(sys.nstates):
        if sys.statePhase[s] != 1:
            continue
        ist = sys.stateStation[s]
        r = sys.stateClass[s]
        if sys.sched[ist] == SchedStrategy.EXT:
            if sys.form == 'J':
                x0[s] = 1.0  # unit mass conserved at the source
            continue
        x0[s] = _class_mass_at_station(sn, ist, r, K)
    sys.x0 = x0


def export_odes_latex(sys, model_name, notation, hide_immediate=False):
    """Render the LaTeX document for the given system.

    Args:
        sys: SymODEs built by solver_fluid_symodes
        model_name: model name shown in the document
        notation: 'scalar' or 'matrix'
        hide_immediate: solver option flag, adds a remark when set

    Returns:
        LaTeX source (str)
    """
    if notation not in ('scalar', 'matrix'):
        raise ValueError("Unknown notation '%s'. Valid notations: scalar, matrix." % notation)
    n = sys.nstates
    L = []
    L.append('% Mean-field fluid ODE system exported by LINE SolverFLD')
    L.append('%% model: %s' % model_name)
    L.append('%% method: %s' % sys.method)
    if sys.form == 'W':
        L.append('% form: dx/dt = W^T*theta(x) + lambda')
    else:
        L.append('% form: dx/dt = J*r(x)')
    L.append('%% notation: %s' % notation)
    L.append('%% nstates: %d' % n)
    if sys.form == 'J':
        L.append('%% nevents: %d' % sys.nevents)
    for s in range(n):
        L.append('%% STATE %d station=%s class=%s phase=%d' % (
            s + 1, sys.stationNames[sys.stateStation[s]], sys.classNames[sys.stateClass[s]], sys.statePhase[s]))
    if sys.form == 'J':
        for e in range(sys.nevents):
            L.append('%% EVENT %d var=%d type=%s coeff=%s' % (
                e + 1, sys.eventVar[e] + 1, sys.factorType[e], _cformat(sys.coeff[e], 15)))

    L.append('\\documentclass{article}')
    L.append('\\usepackage{amsmath}')
    L.append('\\usepackage[margin=2.5cm]{geometry}')
    L.append('\\allowdisplaybreaks')
    L.append('\\setcounter{MaxMatrixCols}{500}')
    L.append('\\begin{document}')
    L.append('\\section*{Mean-field fluid ODE system}')
    L.append('\\noindent Model: \\texttt{%s}. Solver: \\texttt{SolverFLD}, method \\texttt{%s}, %s notation.' % (
        _texesc(model_name), _texesc(sys.method), notation))
    if sys.form == 'W':
        L.append('The system has %d state variables and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}$.' % n)
    else:
        L.append('The system has %d state variables and %d events and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})$.' % (n, sys.nevents))

    L.append('\\subsection*{State variables}')
    L.append('Each state variable $x_{s}$ is the mean number of jobs of a class in a service phase at a station:')
    L.append('\\begin{center}')
    chunk = 48
    for s0 in range(0, n, chunk):
        s1 = min(n, s0 + chunk)
        L.append('\\begin{tabular}{rlll}')
        L.append('\\hline')
        L.append('$s$ & station & class & phase\\\\')
        L.append('\\hline')
        for s in range(s0, s1):
            L.append('%d & \\texttt{%s} & \\texttt{%s} & %d\\\\' % (
                s + 1, _texesc(sys.stationNames[sys.stateStation[s]]),
                _texesc(sys.classNames[sys.stateClass[s]]), sys.statePhase[s]))
        L.append('\\hline')
        L.append('\\end{tabular}')
        if s1 < n:
            L.append('\\par\\medskip')
    L.append('\\end{center}')

    used_stations = sorted(set(int(s) for s in sys.stateStation))
    L.append('\\begin{center}')
    L.append('\\begin{tabular}{rlll}')
    L.append('\\hline')
    L.append('$i$ & station & scheduling & $S_{i}$\\\\')
    L.append('\\hline')
    for i in used_stations:
        L.append('%d & \\texttt{%s} & %s & $%s$\\\\' % (
            i + 1, _texesc(sys.stationNames[i]), _texesc(sys.schedNames[i]), _fmtnum(sys.S[i])))
    L.append('\\hline')
    L.append('\\end{tabular}')
    L.append('\\end{center}')

    T, var_type, var_station, var_class, var_others, const_term = _build_terms(sys)

    defs = _build_defs(sys, var_type, var_station, var_class)
    if defs:
        L.append('\\subsection*{Definitions}')
        L.append('\\begin{align*}')
        L.extend(defs)
        L.append('\\end{align*}')

    if notation == 'scalar':
        L.append('\\subsection*{ODE system (scalar notation)}')
        L.append('\\begin{align}')
        for s in range(n):
            L.append(_render_equation(sys, s, T, var_type, var_station, var_class, var_others,
                                      const_term, s == n - 1))
        L.append('\\end{align}')
    else:
        L.append('\\subsection*{ODE system (matrix notation)}')
        if sys.form == 'W':
            have_lambda = bool(np.any(sys.Alambda != 0))
            L.append('\\begin{equation}')
            if have_lambda:
                L.append('\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}')
            else:
                L.append('\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\,\\theta(\\mathbf{x})')
            L.append('\\end{equation}')
            L.append('with $\\theta_{s}(\\mathbf{x})$ given componentwise by')
            L.append('\\begin{equation*}')
            rows = []
            for v in range(n):
                if sys.isSource[v]:
                    rows.append('0')
                else:
                    rows.append(_factor_tex(v, var_type[v], var_station[v], var_class[v], var_others[v]))
            L.append('\\theta(\\mathbf{x}) = \\begin{bmatrix} %s \\end{bmatrix}' % ' \\\\ '.join(rows))
            L.append('\\end{equation*}')
            L.append('and')
            L.append('\\begin{equation*}')
            L.append('W^{\\top} = %s' % _render_num_matrix(sys.W.T))
            L.append('\\end{equation*}')
            if have_lambda:
                L.append('\\begin{equation*}')
                L.append('\\boldsymbol{\\lambda} = %s^{\\top}' % _render_num_vector(sys.Alambda))
                L.append('\\end{equation*}')
        else:
            L.append('\\begin{equation}')
            L.append('\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})')
            L.append('\\end{equation}')
            L.append('with stoichiometry matrix')
            L.append('\\begin{equation*}')
            J = np.zeros((n, sys.nevents))
            for e in range(sys.nevents):
                J[sys.eventFrom[e], e] -= 1
                J[sys.eventTo[e], e] += 1
            L.append('J = %s' % _render_num_matrix(J))
            L.append('\\end{equation*}')
            L.append('and event rate functions')
            L.append('\\begin{align*}')
            for e in range(sys.nevents):
                fstr = _factor_tex(sys.eventVar[e], sys.factorType[e], sys.factorStation[e],
                                   sys.factorClass[e], sys.factorOthers[e])
                L.append('r_{%d}(\\mathbf{x}) &= %s%s' % (
                    e + 1, _term_tex(sys.coeff[e], fstr), '\\\\' if e < sys.nevents - 1 else ''))
            L.append('\\end{align*}')

    if sys.x0 is not None:
        L.append('\\subsection*{Initial condition}')
        L.append('\\begin{equation*}')
        L.append('\\mathbf{x}(0) = %s^{\\top}' % _render_num_vector(sys.x0))
        L.append('\\end{equation*}')

    L.append('\\subsection*{Remarks}')
    L.append('\\begin{itemize}')
    L.append('\\item For each station $i$, $n_{i}(\\mathbf{x})$ denotes the total mass at the station and $S_{i}$ the number of servers (infinite-server stations use the closed job population, $\\infty$ denotes infinity).')
    L.append('\\item The numerical solver regularizes vanishing denominators with a small positive constant; these regularizations are omitted here.')
    if sys.form == 'J':
        if any(t in ('fcfsw', 'fcfsws') for t in sys.factorType):
            L.append('\\item At FCFS stations, the mean phase residence times $w_{u} = -1/[D_{0}]_{kk}$ weight the backlog $\\hat{n}_{i}$; the factors $w_{u}$ of the departing phases are folded into the rate coefficients.')
        if any(t == 'dpsmin' for t in sys.factorType):
            L.append('\\item At DPS stations, weights are normalized to sum to one and the weight $w_{ir}$ of the departing class is folded into the rate coefficient; the class shares $w_{ir}x/\\tilde{n}_{i}$ divide the station capacity $\\min(n_{i},S_{i})$, so they sum to one whenever the station is busy.')
    any_fcfs = any(sys.sched[sys.stateStation[s]] == SchedStrategy.FCFS for s in range(n))
    if any_fcfs and sys.method in ('matrix', 'closing'):
        L.append('\\item For FCFS stations with non-exponential service, the solver may iteratively re-fit the service distributions (non-exponential approximation); the exported system uses the nominal model parameters.')
    if hide_immediate:
        L.append('\\item \\texttt{hide\\_immediate} is enabled in the solver options: the numerical integration may further eliminate immediate transitions by state-space reduction; the exported system is the unreduced one.')
    L.append('\\end{itemize}')
    L.append('\\end{document}')

    return '\n'.join(L) + '\n'


def _build_terms(sys):
    n = sys.nstates
    T = np.zeros((n, n))
    var_type = [None] * n
    var_station = [0] * n
    var_class = [0] * n
    var_others = [None] * n
    const_term = np.zeros(n)
    if sys.form == 'W':
        T = sys.W.T.copy()
        T[:, sys.isSource] = 0.0  # theta of Source states is identically zero
        for v in range(n):
            if not sys.isSource[v]:
                i = sys.stateStation[v]
                var_type[v] = 'lin' if np.isinf(sys.S[i]) else sys.smoothing
                var_station[v] = i
                var_class[v] = sys.stateClass[v]
        const_term = sys.Alambda.copy()
    else:
        for e in range(sys.nevents):
            v = sys.eventVar[e]
            T[sys.eventFrom[e], v] -= sys.coeff[e]
            T[sys.eventTo[e], v] += sys.coeff[e]
            if var_type[v] is None:
                var_type[v] = sys.factorType[e]
                var_station[v] = sys.factorStation[e]
                var_class[v] = sys.factorClass[e]
                var_others[v] = sys.factorOthers[e]
    return T, var_type, var_station, var_class, var_others, const_term


def _build_defs(sys, var_type, var_station, var_class):
    defs = []
    n = sys.nstates
    M = len(sys.stationNames)
    K = len(sys.classNames)
    need_n = [False] * M
    need_nt = [False] * M
    need_nh = [False] * M
    gdef = [None] * M
    for v in range(n):
        f = var_type[v]
        if f is None:
            continue
        i = var_station[v]
        if f == 'min':
            need_n[i] = True
            gdef[i] = 'g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{n_{%d}(\\mathbf{x})}' % (
                i + 1, i + 1, _fmtnum(sys.S[i]), i + 1)
        elif f == 'pnorm':
            need_n[i] = True
            gdef[i] = 'g_{%d}(\\mathbf{x}) &= \\Bigl(1 + \\bigl(n_{%d}(\\mathbf{x})/%s\\bigr)^{%s}\\Bigr)^{-1/%s}' % (
                i + 1, i + 1, _fmtnum(sys.S[i]), _fmtnum(sys.pstar[i]), _fmtnum(sys.pstar[i]))
        elif f == 'dpsmin':
            need_n[i] = True
            need_nt[i] = True
            gdef[i] = 'g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\tilde{n}_{%d}(\\mathbf{x})}' % (
                i + 1, i + 1, _fmtnum(sys.S[i]), i + 1)
        elif f == 'dpspw':
            need_n[i] = True
            need_nt[i] = True
        elif f in ('fcfsw', 'fcfsws'):
            need_n[i] = True
            need_nh[i] = True
            if f == 'fcfsw':
                gdef[i] = 'g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\hat{n}_{%d}(\\mathbf{x})}' % (
                    i + 1, i + 1, _fmtnum(sys.S[i]), i + 1)
            else:
                gdef[i] = 'g_{%d}(\\mathbf{x}) &= \\frac{\\mathrm{softmin}\\bigl(n_{%d}(\\mathbf{x}),\\, %s\\bigr)}{\\hat{n}_{%d}(\\mathbf{x})}' % (
                    i + 1, i + 1, _fmtnum(sys.S[i]), i + 1)
    for i in range(M):
        if need_n[i]:
            terms = ' + '.join('x_{%d}' % (v + 1) for v in range(n) if sys.stateStation[v] == i)
            defs.append('n_{%d}(\\mathbf{x}) &= %s\\\\' % (i + 1, terms))
    for i in range(M):
        if need_nt[i]:
            parts = []
            for r in range(K):
                terms = ' + '.join('x_{%d}' % (v + 1) for v in range(n)
                                   if sys.stateStation[v] == i and sys.stateClass[v] == r)
                if terms:
                    parts.append('%s\\,(%s)' % (_fmtnum(sys.dpsw[i, r]), terms))
            defs.append('\\tilde{n}_{%d}(\\mathbf{x}) &= %s\\\\' % (i + 1, ' + '.join(parts)))
    for i in range(M):
        if need_nh[i]:
            parts = ' + '.join('%s\\,x_{%d}' % (_fmtnum(sys.fcfsPhaseW[v]), v + 1)
                               for v in range(n) if sys.stateStation[v] == i)
            defs.append('\\hat{n}_{%d}(\\mathbf{x}) &= %s\\\\' % (i + 1, parts))
    for i in range(M):
        if gdef[i] is not None:
            defs.append(gdef[i] + '\\\\')
    for v in range(n):
        if var_type[v] == 'dpspw':
            i = var_station[v]
            r = var_class[v]
            d = ('g_{%d,%d}(\\mathbf{x}) &= \\begin{cases} 1 & n_{%d}(\\mathbf{x}) \\le %s\\\\ '
                 '\\dfrac{%s}{\\tilde{n}_{%d}(\\mathbf{x})} & n_{%d}(\\mathbf{x}) > %s \\end{cases}\\\\') % (
                i + 1, r + 1, i + 1, _fmtnum(sys.S[i]), _fmtnum(sys.S[i] * sys.dpsw[i, r]),
                i + 1, i + 1, _fmtnum(sys.S[i]))
            if d not in defs:
                defs.append(d)
    if any(var_type[v] == 'fcfsws' for v in range(n)):
        defs.append('\\mathrm{softmin}(a,b) &= \\frac{a\\,e^{-\\alpha a} + b\\,e^{-\\alpha b}}{e^{-\\alpha a} + e^{-\\alpha b}}, \\qquad \\alpha = %s\\\\' % _fmtnum(sys.alpha))
    if defs and defs[-1].endswith('\\\\'):
        defs[-1] = defs[-1][:-2]
    return defs


def _render_equation(sys, sidx, T, var_type, var_station, var_class, var_others, const_term, is_last):
    terms = []
    for v in range(sys.nstates):
        c = T[sidx, v]
        if c != 0:
            terms.append((c, _factor_tex(v, var_type[v], var_station[v], var_class[v], var_others[v])))
    if const_term[sidx] != 0:
        terms.append((const_term[sidx], ''))
    if not terms:
        rhs = '0'
    else:
        parts = []
        for k, (c, fstr) in enumerate(terms):
            body = _term_tex(abs(c), fstr)
            if k == 0:
                parts.append(('-' if c < 0 else '') + body)
            else:
                parts.append((' - ' if c < 0 else ' + ') + body)
            if (k + 1) % 4 == 0 and k < len(terms) - 1:
                parts.append('\\nonumber\\\\\n&\\quad ')
        rhs = ''.join(parts)
    return '\\frac{\\mathrm{d}x_{%d}}{\\mathrm{d}t} &= %s%s' % (sidx + 1, rhs, '' if is_last else '\\\\')


def _factor_tex(v, ftype, station, class_idx, others):
    if ftype == 'lin':
        return 'x_{%d}' % (v + 1)
    if ftype in ('min', 'pnorm', 'dpsmin', 'fcfsw', 'fcfsws'):
        return 'x_{%d}\\,g_{%d}(\\mathbf{x})' % (v + 1, station + 1)
    if ftype == 'dpspw':
        return 'x_{%d}\\,g_{%d,%d}(\\mathbf{x})' % (v + 1, station + 1, class_idx + 1)
    # ext1
    if not others:
        return ''  # single-phase source class: constant unit mass
    return '\\bigl(1 - %s\\bigr)' % ' - '.join('x_{%d}' % (u + 1) for u in others)


def _term_tex(c, fstr):
    if not fstr:
        return _fmtnum(c)
    if c == 1:
        return fstr
    return '%s\\,%s' % (_fmtnum(c), fstr)


def _render_num_matrix(A):
    A = np.atleast_2d(np.asarray(A, dtype=float))
    rows = [' & '.join(_fmtnum(A[r, c]) for c in range(A.shape[1])) for r in range(A.shape[0])]
    body = ' \\\\ '.join(rows)
    if max(A.shape) > 12:
        return '{\\scriptsize\\begin{bmatrix} %s \\end{bmatrix}}' % body
    return '\\begin{bmatrix} %s \\end{bmatrix}' % body


def _render_num_vector(v):
    body = ' & '.join(_fmtnum(x) for x in np.asarray(v, dtype=float).flatten())
    return '\\begin{pmatrix} %s \\end{pmatrix}' % body


def _fmtnum(v):
    """Compact LaTeX-safe number formatting, matching MATLAB's %.8g style."""
    v = float(v)
    if np.isinf(v):
        return '\\infty' if v > 0 else '-\\infty'
    if v == round(v) and abs(v) < 1e15:
        return '%d' % int(v)
    return '%.8g' % v


def _cformat(v, prec):
    """C-style %g formatting used for parity with MATLAB sprintf."""
    v = float(v)
    if v == round(v) and abs(v) < 1e15:
        return '%d' % int(v)
    return ('%.' + str(prec) + 'g') % v


def _texesc(s):
    import re
    return re.sub(r'([_%&#])', r'\\\1', str(s))


# Symbolic drift, for the computer algebra backend.
# see _kb/06-solver-catalog.md (Fluid: "Symbolic drift and Jacobian") for which
# factor types are smooth and exportable vs which carry a min/branch.
_SMOOTH_FACTORS = ('lin', 'ext1', 'fcfsws')


def state_variables(sys):
    """State variable names of the exported drift, x1 ... xn.

    Args:
        sys: SymODEs instance

    Returns:
        list of variable names
    """
    return ['x%d' % (s + 1) for s in range(sys.nstates)]


def symbolic_drift(sys):
    """Right-hand side of the ODE system as expression strings, one per state
    variable, in the format the symbolic backend parses.

    ONLY SMOOTH DRIFTS ARE EXPORTED. The default, matrix, closing and statedep
    methods scale rates by min(n_i, S_i), which is not differentiable at
    n_i = S_i, so their Jacobian does not exist there; emitting a one-sided
    derivative would be a silent lie exactly at the regime switch that
    matters. Use the p-norm smoothing (options.pstar, method matrix or pnorm)
    or the softmin method.

    FineTol is carried in exactly the places the integrated systems put it,
    and nowhere else: the p-norm drift offsets the station total, the softmin
    drift offsets the phase-weighted total but not the plain station total
    that feeds the softmin, and the closing rates offset neither. Mirrors
    MATLAB @SolverFLD/getSymbolicDrift.m and the JAR FluidODEsExporter.

    Args:
        sys: SymODEs instance

    Returns:
        list of expression strings, one per state variable

    Raises:
        ValueError: if the drift is not differentiable
    """
    from ...constants import GlobalConstants

    variables = state_variables(sys)
    eps0 = _num(GlobalConstants.FineTol)
    if sys.form == 'W':
        if sys.smoothing != 'pnorm':
            raise ValueError(
                'The drift of this method scales rates by min(n_i, S_i), which is not '
                'differentiable at n_i = S_i, so it has no Jacobian there. Set '
                "options.pstar to use the p-norm smoothing, or use the 'softmin' "
                'method.')
        return _wform_drift(sys, variables, eps0)
    if sys.form == 'J':
        return _jform_drift(sys, variables, eps0)
    raise ValueError("unsupported ODE form '%s'" % sys.form)


def _wform_drift(sys, variables, eps0):
    """dx/dt = W' * theta(x) + Alambda, with the p-norm smoothed
        theta_s = x_s / (1 + (n_i/S_i)^p_i)^(1/p_i),  theta_s = 0 at a Source,
    mirroring pnorm_ode in solver_fluid_matrix.
    """
    n = sys.nstates
    theta = [None] * n
    for s in range(n):
        if sys.isSource[s]:
            theta[s] = '0'
            continue
        i = int(sys.stateStation[s])
        ni = _station_sum(sys, i, variables, eps0)
        S = float(sys.S[i])
        p = float(sys.pstar[i])
        if S <= 0 or p <= 0:
            theta[s] = variables[s]
        else:
            theta[s] = '%s/(1 + (%s/%s)^%s)^(1/%s)' % (
                variables[s], ni, _num(S), _num(p), _num(p))

    rhs = []
    for s in range(n):
        terms = []
        for t in range(n):
            w = float(sys.W[t, s])
            if w == 0 or theta[t] == '0':
                continue
            terms.append('(%s)*(%s)' % (_num(w), theta[t]))
        if sys.Alambda[s] != 0:
            terms.append(_num(sys.Alambda[s]))
        rhs.append(' + '.join(terms) if terms else '0')
    return rhs


def _jform_drift(sys, variables, eps0):
    """dx/dt = J * r(x), with r_e = coeff(e) * factor_e(x). Only the smooth
    factor types are exportable: 'min' (PS/FCFS under closing and statedep),
    'fcfsw' (statedep FCFS), 'dpsmin' (closing DPS) and 'dpspw' (piecewise DPS)
    all carry a min or a branch.
    """
    rate = [None] * sys.nevents
    for e in range(sys.nevents):
        ftype = sys.factorType[e]
        if ftype not in _SMOOTH_FACTORS:
            raise ValueError(
                "Event %d scales its rate by the non-smooth factor '%s', which has no "
                'derivative where the regime switches, so the system has no Jacobian. '
                "Use the 'softmin' method, or the p-norm smoothing of the 'matrix' "
                'method.' % (e + 1, ftype))
        v = variables[int(sys.eventVar[e])]
        if ftype == 'lin':
            factor = v
        elif ftype == 'ext1':
            # 1 - sum of the class's phases 2..end at the source
            others = sys.factorOthers[e]
            if not others:
                factor = '1'
            else:
                factor = '(1 - (%s))' % ' + '.join(variables[int(k)] for k in others)
        else:  # fcfsws
            i = int(sys.factorStation[e])
            # ode_softmin: ni is the raw station total, wni carries FineTol.
            ni = _station_sum(sys, i, variables, '0')
            nhat = _phase_weighted_station_sum(sys, i, variables, eps0)
            factor = '%s*(%s)/(%s)' % (v, _softmin_expr(ni, _num(sys.S[i]), sys.alpha), nhat)
        rate[e] = '(%s)*(%s)' % (_num(sys.coeff[e]), factor)

    rhs = []
    for s in range(sys.nstates):
        terms = []
        for e in range(sys.nevents):
            j = 0.0
            if sys.eventFrom[e] == s:
                j -= 1.0
            if sys.eventTo[e] == s:
                j += 1.0
            if j == 0:
                continue
            terms.append('(%s)*(%s)' % (_num(j), rate[e]))
        rhs.append(' + '.join(terms) if terms else '0')
    return rhs


def _softmin_expr(x, y, alpha):
    """Smooth minimum in its weighted-average form,
      (x e^{-a x} + y e^{-a y}) / (e^{-a x} + e^{-a y}),
    which is what softmin computes; the implementation rewrites it as
    lo + gap*w/(1+w) only to keep the exponent argument non-positive, an
    overflow guard that is meaningless symbolically and would introduce the
    min/max branch this export exists to avoid.
    """
    a = _num(alpha)
    return ('((%s)*exp(-(%s)*(%s)) + (%s)*exp(-(%s)*(%s)))'
            '/(exp(-(%s)*(%s)) + exp(-(%s)*(%s)))'
            % (x, a, x, y, a, y, a, x, a, y))


def _station_sum(sys, i, variables, offset):
    """Total fluid mass at station i, plus the offset its consumer uses."""
    parts = [variables[k] for k in range(sys.nstates) if sys.stateStation[k] == i]
    return _join_sum(parts, offset)


def _phase_weighted_station_sum(sys, i, variables, offset):
    """sum_u w_u x_u over the states of station i, with w_u the mean phase
    time weights (nhat in ode_softmin)."""
    parts = []
    for k in range(sys.nstates):
        if sys.stateStation[k] == i:
            wt = float(sys.fcfsPhaseW[k])
            if wt != 0:
                parts.append('(%s)*%s' % (_num(wt), variables[k]))
    return _join_sum(parts, offset)


def _join_sum(parts, offset):
    """Sum of PARTS, with OFFSET added only when it is not the literal zero."""
    if not parts:
        return '(%s)' % offset
    if offset == '0':
        return '(%s)' % ' + '.join(parts)
    return '(%s + %s)' % (offset, ' + '.join(parts))


def _num(v):
    """Decimal text the symbolic backend reads as an exact rational."""
    v = float(v)
    if v == round(v) and abs(v) < 1e15:
        return '%d' % int(round(v))
    return '%.17g' % v
