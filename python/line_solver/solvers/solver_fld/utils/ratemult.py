"""
Time-varying rate multipliers for the fluid ODE.

Port of MATLAB fluid_interpcols.m and solver_fluid_ratemult.m.

The fluid closing ODE is autonomous: every event rate is a fixed base rate
scaled by a state-dependent factor. A time-varying input (a non-homogeneous
Poisson source, or a caller-supplied demand trajectory) is expressed as a
per-event multiplicative factor m(t), because the closing rate is
rate = rateBase .* theta(x) and rateBase is linear in the station-class
service/arrival rate. Evaluating m(t) by clamped piecewise-linear
interpolation (fluid_interpcols) makes the otherwise autonomous ODE
non-autonomous while leaving the legacy form untouched when no channel is
configured.

Three independent, composable channels are honoured, all read from
options.config:

  rate_traj   (tgrid, Mmat) - a caller-supplied event multiplier matrix
              (numEvents x len(tgrid)); used by the coupled LN layer transient.
  nhpp_sched  a sequence of entries with fields station, class, nhpp, where
              nhpp exposes getRateAt/getBreakpoints/getPeriod/isCyclic. The
              nominal (time-average) rate baked into the base rate is the
              station-class rate, and the multiplier is getRateAt(t)/nominal.
  rate_sched  a sequence of entries with fields station, class, tgrid, rates
              and optionally nominal; an explicit per-(station,class) rate
              trajectory.
"""

import numpy as np


def fluid_interpcols(tg, B, tt):
    """Clamped piecewise-linear interpolation of the columns of B at time tt.

    B is (nrows x ngrid) with column j sampled at time tg[j]; tg is strictly
    increasing. Times outside [tg[0], tg[-1]] are clamped to the boundary
    columns (zero-order hold outside the grid). This is the shared
    time-varying-input evaluator used both by the TBI transient (frozen
    complement trajectory) and by the rate-multiplier closure.

    Args:
        tg: (ngrid,) strictly increasing time grid.
        B: (nrows x ngrid) column samples.
        tt: scalar time.

    Returns:
        (nrows,) interpolated column.
    """
    tg = np.asarray(tg, dtype=float).ravel()
    B = np.asarray(B, dtype=float)
    if B.ndim == 1:
        B = B.reshape(1, -1)
    if tt <= tg[0]:
        return B[:, 0].copy()
    if tt >= tg[-1]:
        return B[:, -1].copy()
    j = int(np.searchsorted(tg, tt, side='right')) - 1
    wj = (tt - tg[j]) / (tg[j + 1] - tg[j])
    return (1.0 - wj) * B[:, j] + wj * B[:, j + 1]


def _field(entry, name, default=None):
    """Read a field of a schedule entry given as a dict or as an object."""
    if isinstance(entry, dict):
        if name in entry:
            return entry[name]
        if name == 'class' and 'jobclass' in entry:
            return entry['jobclass']
        return default
    if hasattr(entry, name):
        return getattr(entry, name)
    if name == 'class' and hasattr(entry, 'jobclass'):
        return getattr(entry, 'jobclass')
    return default


def _config(options):
    """The options.config bag as a dict (empty when absent)."""
    cfg = getattr(options, 'config', None)
    if cfg is None:
        return {}
    if isinstance(cfg, dict):
        return cfg
    return {k: getattr(cfg, k) for k in dir(cfg) if not k.startswith('_')}


def _horizon(options):
    """Integration horizon (t0, tend); tend is inf when unbounded."""
    t0 = 0.0
    tend = np.inf
    ts = getattr(options, 'timespan', None)
    if ts is not None and len(ts) >= 2:
        if np.isfinite(ts[0]):
            t0 = float(ts[0])
        tend = float(ts[1])
    return t0, tend


def nhpp_steps(nh, t0, thi):
    """Step-faithful (time, rate) sampling of a piecewise-constant intensity.

    Each segment of the schedule contributes two samples, at its start and just
    before its end, so clamped-linear interpolation reproduces the step with a
    negligible transition ramp. A cyclic schedule is unrolled over [t0, thi].

    Args:
        nh: process exposing getBreakpoints/getPeriod/isCyclic/getRateAt.
        t0: horizon start.
        thi: horizon end.

    Returns:
        (seg_t, seg_r) strictly increasing times and the rate at each.
    """
    bp = np.asarray(nh.getBreakpoints(), dtype=float).ravel()
    period = float(nh.getPeriod())
    if nh.isCyclic() and np.isfinite(period) and period > 0:
        kmax = int(np.ceil((thi - t0) / period)) + 2
        bounds = np.concatenate([bp + k * period for k in range(-1, kmax + 1)])
    else:
        bounds = bp
    inner = bounds[(bounds > t0) & (bounds < thi)]
    bounds = np.unique(np.concatenate([[t0], inner, [thi]]))
    neps = max(1e-9, 1e-6 * (thi - t0))
    nb = bounds.size - 1
    seg_t = np.zeros(2 * nb)
    seg_r = np.zeros(2 * nb)
    for k in range(nb):
        a = bounds[k]
        b = bounds[k + 1]
        r = float(nh.getRateAt(0.5 * (a + b)))
        seg_t[2 * k] = a
        seg_r[2 * k] = r
        seg_t[2 * k + 1] = max(a + neps, b - neps)
        seg_r[2 * k + 1] = r
    return seg_t, seg_r


def merge_multipliers(tg1, M1, tg2, M2, nrows):
    """Merge two multiplier trajectories on the union grid by elementwise product.

    Identity where a source is silent, so a missing channel leaves the other
    unchanged.

    Args:
        tg1, M1: first trajectory (either may be None/empty).
        tg2, M2: second trajectory (either may be None/empty).
        nrows: number of rows of the merged matrix.

    Returns:
        (tgrid, Mmat), both None when both inputs are empty.
    """
    empty1 = M1 is None or np.size(M1) == 0
    empty2 = M2 is None or np.size(M2) == 0
    if empty1:
        return (None, None) if empty2 else (tg2, M2)
    if empty2:
        return tg1, M1
    tg = np.unique(np.concatenate([np.asarray(tg1, dtype=float).ravel(),
                                   np.asarray(tg2, dtype=float).ravel()]))
    Mg = np.ones((nrows, tg.size))
    for j in range(tg.size):
        Mg[:, j] = fluid_interpcols(tg1, M1, tg[j]) * fluid_interpcols(tg2, M2, tg[j])
    return tg, Mg


def ratemult_entries(enabled, nominal, options):
    """Collect the per-(station,class) rate multiplier trajectories.

    Reads the nhpp_sched and rate_sched channels of options.config. The
    rate_traj channel is event-indexed and is handled by
    solver_fluid_ratemult instead.

    Args:
        enabled: (M x K) bool, whether a class is served at a station.
        nominal: (M x K) float, the nominal rate baked into the base rate.
        options: solver options carrying config and timespan.

    Returns:
        dict with keys 'nhpp' and 'sched', each a list of
        (station, class, seg_t, rowmult) tuples.
    """
    cfg = _config(options)
    enabled = np.asarray(enabled, dtype=bool)
    nominal = np.asarray(nominal, dtype=float)
    out = {'nhpp': [], 'sched': []}

    sched = cfg.get('nhpp_sched')
    if sched is not None and len(sched) > 0:
        t0, tend = _horizon(options)
        for entry in sched:
            i = int(_field(entry, 'station'))
            c = int(_field(entry, 'class'))
            nh = _field(entry, 'nhpp')
            if nh is None or not enabled[i, c]:
                continue
            nom = float(nominal[i, c])
            if not (nom > 0):
                continue
            # horizon: for a non-finite timespan use a few periods so a cyclic
            # schedule is represented rather than clamped after one segment
            period = float(nh.getPeriod())
            if not np.isfinite(tend):
                if np.isfinite(period) and period > 0:
                    thi = t0 + 3.0 * period
                else:
                    thi = t0 + 1.0
            else:
                thi = tend
            seg_t, seg_r = nhpp_steps(nh, t0, thi)
            out['nhpp'].append((i, c, seg_t, seg_r / nom))

    rs = cfg.get('rate_sched')
    if rs is not None and len(rs) > 0:
        for entry in rs:
            i = int(_field(entry, 'station'))
            c = int(_field(entry, 'class'))
            if not enabled[i, c]:
                continue
            nom = _field(entry, 'nominal')
            if nom is None or np.size(nom) == 0:
                nom = nominal[i, c]
            nom = float(np.asarray(nom).ravel()[0])
            if not (nom > 0):
                continue
            seg_t = np.asarray(_field(entry, 'tgrid'), dtype=float).ravel()
            seg_r = np.asarray(_field(entry, 'rates'), dtype=float).ravel()
            out['sched'].append((i, c, seg_t, seg_r / nom))

    return out


def expand_entries_rows(entries, nrows, rows_of):
    """Merge per-(station,class) entries into one (nrows x ngrid) trajectory.

    Args:
        entries: list of (station, class, seg_t, rowmult).
        nrows: number of rows of the result.
        rows_of: callable (station, class) -> boolean row mask.

    Returns:
        (tgrid, Mmat), both None when there is no entry.
    """
    tg = None
    Mm = None
    for (i, c, seg_t, rowmult) in entries:
        rows = rows_of(i, c)
        if not np.any(rows):
            continue
        thisM = np.ones((nrows, seg_t.size))
        thisM[rows, :] = np.tile(rowmult, (int(np.count_nonzero(rows)), 1))
        tg, Mm = merge_multipliers(tg, Mm, seg_t, thisM, nrows)
    return tg, Mm


def solver_fluid_ratemult(num_events, enabled, q_indices, Kic, nominal,
                          event_idx, options):
    """Build the time-varying per-event rate multiplier for the closing ODE.

    Port of MATLAB solver_fluid_ratemult.m. Event e is scaled at time t by
    fluid_interpcols(tgrid, Mmat, t)[e].

    Args:
        num_events: number of events in the closing ODE.
        enabled: (M x K) bool, whether a class is served at a station.
        q_indices: (M x K) starting state index of each station-class.
        Kic: (M x K) number of phases of each station-class.
        nominal: (M x K) nominal rate baked into the base rate.
        event_idx: (num_events,) state index driving each event.
        options: solver options carrying config and timespan.

    Returns:
        (tgrid, Mmat) with Mmat of shape (num_events x len(tgrid)), or
        (None, None) when no time-varying channel is configured, so the caller
        keeps the legacy autonomous closure.
    """
    cfg = _config(options)
    if not cfg:
        return None, None

    event_idx = np.asarray(event_idx).ravel()
    q_indices = np.asarray(q_indices)
    Kic = np.asarray(Kic)

    def rows_of(i, c):
        rows = np.zeros(num_events, dtype=bool)
        for kic in range(int(Kic[i, c])):
            rows |= (event_idx == (q_indices[i, c] + kic))
        return rows

    # source (1): caller-supplied rate_traj
    user_tg = None
    user_M = None
    traj = cfg.get('rate_traj')
    if traj is not None and len(traj) >= 2:
        user_tg = np.asarray(traj[0], dtype=float).ravel()
        user_M = np.asarray(traj[1], dtype=float)
        if user_M.shape[0] != num_events:
            raise ValueError(
                'rate_traj multiplier matrix has %d rows but the closing ODE '
                'has %d events.' % (user_M.shape[0], num_events))

    entries = ratemult_entries(enabled, nominal, options)

    # source (2): NHPP source intensities
    nhpp_tg, nhpp_M = expand_entries_rows(entries['nhpp'], num_events, rows_of)
    # source (3): explicit per-(station,class) rate trajectories
    sched_tg, sched_M = expand_entries_rows(entries['sched'], num_events, rows_of)

    tgrid, Mmat = merge_multipliers(user_tg, user_M, nhpp_tg, nhpp_M, num_events)
    tgrid, Mmat = merge_multipliers(tgrid, Mmat, sched_tg, sched_M, num_events)
    return tgrid, Mmat


def ratemult_max_step(tgrid):
    """Integration step cap that resolves the steps of a multiplier grid.

    A piecewise-constant schedule is a sequence of near-discontinuities; an
    adaptive integrator with an unbounded step can walk straight over a short
    segment. Cap the step at a quarter of the shortest genuine segment, the
    paired samples that encode one step (separated by the sampling epsilon)
    being ignored.

    Args:
        tgrid: the multiplier time grid, or None.

    Returns:
        The step cap, or inf when there is nothing to resolve.
    """
    if tgrid is None:
        return np.inf
    tgrid = np.asarray(tgrid, dtype=float).ravel()
    if tgrid.size < 2:
        return np.inf
    span = tgrid[-1] - tgrid[0]
    if not (span > 0):
        return np.inf
    dt = np.diff(tgrid)
    dt = dt[dt > 1e-5 * span]
    if dt.size == 0:
        return np.inf
    return float(np.min(dt)) / 4.0


__all__ = [
    'fluid_interpcols',
    'nhpp_steps',
    'merge_multipliers',
    'ratemult_entries',
    'expand_entries_rows',
    'solver_fluid_ratemult',
    'ratemult_max_step',
]
