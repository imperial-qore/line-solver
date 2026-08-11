"""Native SolverBA: bound-analysis solver.

Subclasses the native SolverMVA to reuse its model parsing (rates, demands,
njobs, nservers, sched, chain info), and replaces runAnalyzer with the
bound-analysis dispatch in solver_ba_analyzer. Mirrors matlab/src/solvers/BA.
"""

import numpy as np

from ..solver_mva.solver_mva import SolverMVA as _NativeMVA
from .solver_ba_analyzer import (
    solver_ba_analyzer, BA_ASYM, BA_CHAIN, BA_HIER, BA_LR, BA_QRF_LP,
    BA_QRF_NLP)

# QRF reduction bounds are served under SolverBA in MATLAB/JAR, and now
# natively as well. Both families route through solver_ctmc_qrf_analyzer:
#   qrf.bas, qrf.rsrd            -- LP, via the qr_bounds_* backends (HiGHS).
#   qr, qrf.mmi, qrf.mem,        -- NLP, via api.mapqn.qrf_noblo_*.
#   qrf.mmi.ld, qrf.mmi.linear
# The NLP tokens were withheld until 2026-07-20 pending end-to-end coverage of
# the adapter path. Exercising it found a transposed v in the adapter and a
# throughput inversion valid only for a delay reference station; with both
# fixed, K == 1 reproduces the CTMC exactly and K == 2 lands inside the glpsol
# bound range. Infinite-server stations are rejected: qrf_noblo_* models every
# station as a single server.
BA_METHODS = ['default'] + sorted(
    BA_ASYM | BA_CHAIN | BA_HIER | BA_LR | BA_QRF_LP | BA_QRF_NLP)


class SolverBA(_NativeMVA):
    """Bound-analysis solver (asymptotic + hierarchical throughput/queue-length
    bounds, plus the QRF LP-based bounds).  cub is upper-only; mbjb.lower is the
    multiclass BJB lower bound (Kerola eq. 10)."""

    def __init__(self, model, method_or_options=None, **kwargs):
        level = kwargs.pop('level', None)
        super().__init__(model, method_or_options, **kwargs)
        if level is not None:
            setattr(self.options, 'level', int(level))
        elif not hasattr(self.options, 'level'):
            setattr(self.options, 'level', 2)

    def getName(self):
        return 'SolverBA'

    def list_valid_methods(self):
        return list(BA_METHODS)

    listValidMethods = list_valid_methods

    def runAnalyzer(self):
        method = str(getattr(self, 'method', 'default')).lower()
        self._result = solver_ba_analyzer(self, method, self.options)
        return self

    run_analyzer = runAnalyzer

    def get_bounds(self):
        """Return {Tlower,Tupper,Qlower,Qupper} for the current method's family.
        One-sided families (cub upper-only, mbjb/ldbcmp lower-only) return NaN on
        the missing side.

        Each side is re-run through self.options, so the caller's full option set
        (notably `level`) is inherited: a hierarchical family reached through
        get_bounds must tighten as level is raised, not silently revert to the
        default level of 2."""
        fam = str(getattr(self, 'method', 'default')).lower().split('.')[0]
        valid = self.list_valid_methods()
        out = {'Tlower': np.nan, 'Tupper': np.nan,
               'Qlower': np.nan, 'Qupper': np.nan}
        for side, kt, kq in (('lower', 'Tlower', 'Qlower'),
                             ('upper', 'Tupper', 'Qupper')):
            m = fam + '.' + side
            if m in valid:
                r = solver_ba_analyzer(self, m, self.options)
                out[kt] = r['TN']
                out[kq] = r['QN']
        return out

    getBounds = get_bounds

    def _expand_bound(self, x, M, K):
        """Normalize a bracket side to an (M,K) array. A one-sided family leaves
        its missing side as the scalar NaN returned by get_bounds; expand it so
        the table keeps full shape with NaN entries."""
        a = np.asarray(x, dtype=float)
        if a.size == 0:
            return np.full((M, K), np.nan)
        if a.ndim == 0:
            return np.full((M, K), float(a))
        return a.reshape(M, K) if a.size == M * K else a

    def get_bounds_table(self, keepDisabled=False):
        """Return the {lower,upper} bracket per station and class as a DataFrame,
        in the layout of getAvgTable. Columns:
        Station, JobClass, Qlower, Qupper, Tlower, Tupper.

        One-sided families (cub upper-only, mbjb/ldbcmp lower-only) carry NaN on
        the missing side; NaN is preserved, never replaced by zero."""
        import pandas as pd
        b = self.get_bounds()
        M = self.nstations
        K = self.nclasses
        Ql = self._expand_bound(b['Qlower'], M, K)
        Qu = self._expand_bound(b['Qupper'], M, K)
        Tl = self._expand_bound(b['Tlower'], M, K)
        Tu = self._expand_bound(b['Tupper'], M, K)
        station_names = list(getattr(self, 'station_names', []))
        class_names = list(getattr(self, 'class_names', []))
        rows = []
        for i in range(M):
            sname = station_names[i] if i < len(station_names) else 'Station%d' % i
            for k in range(K):
                cname = class_names[k] if k < len(class_names) else 'Class%d' % k
                vals = np.array([Ql[i, k], Qu[i, k], Tl[i, k], Tu[i, k]], dtype=float)
                present = vals[~np.isnan(vals)]
                # Mirror getAvgTable's drop of disabled station-class pairs, but
                # NaN-safe: keep the row when any value that is present is
                # nonzero, so an all-NaN side never removes the row.
                if keepDisabled or present.size == 0 or np.any(present != 0):
                    rows.append({
                        'Station': sname,
                        'JobClass': cname,
                        'Qlower': Ql[i, k],
                        'Qupper': Qu[i, k],
                        'Tlower': Tl[i, k],
                        'Tupper': Tu[i, k],
                    })
        return pd.DataFrame(
            rows,
            columns=['Station', 'JobClass', 'Qlower', 'Qupper', 'Tlower', 'Tupper'])

    getBoundsTable = get_bounds_table
