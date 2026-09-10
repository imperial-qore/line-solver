"""
Patience (time-to-abandon) handles derived from a NetworkStruct.

Native Python twin of matlab/src/api/sn/sn_patience_handles.m. The abandonment
solvers need the patience law as FUNCTIONS -- a complementary cdf and a hazard
rate -- not as moments, because that is what the underlying theory consumes:
Whitt's engineering solution reads the hazard near the origin, and the fluid
models integrate the ccdf. LINE stores the law as a MAP/PH pair, from which both
are available in closed form.
"""

from typing import Any, Callable, Dict, Optional

import numpy as np


def sn_patience_handles(sn, ist: int, r: int) -> Optional[Dict[str, Any]]:
    """
    Build ccdf, pdf and hazard handles for the patience law of station ``ist``,
    class ``r``.

    Args:
        sn: the NetworkStruct
        ist: station index
        r: class index

    Returns:
        Dict with ``ccdf``, ``pdf``, ``hazard`` (callables), ``mean``,
        ``isExponential`` and ``rate``; ``None`` when the station-class pair has
        no reneging patience configured.

    See also:
        matlab/src/api/sn/sn_patience_handles.m
    """
    from ...lang.base import ImpatienceType
    from ...constants import ProcessType

    cls = getattr(sn, 'impatienceClass', None)
    if cls is None or int(np.asarray(cls)[ist, r]) != int(ImpatienceType.RENEGING):
        return None
    proc = getattr(sn, 'patienceProc', None) or getattr(sn, 'impatienceProc', None)
    pair = None
    if proc is not None:
        try:
            pair = proc[ist][r]
        except Exception:
            pair = None
    rate = float(np.asarray(getattr(sn, 'impatienceMu', np.zeros((1, 1))))[ist, r])
    ptype = int(np.asarray(getattr(sn, 'impatienceType', np.zeros((1, 1), dtype=int)))[ist, r])
    # ProcessType is an Enum here, an integer id in the struct: compare on
    # the value, since the enum numbering differs from MATLAB's anyway.
    is_exp = (ptype == ProcessType.EXP.value)

    if pair is None:
        if rate <= 0:
            return None
        # Only the rate is on record, so the law is exponential by construction.
        return {
            'ccdf': (lambda t, _r=rate: float(np.exp(-_r * np.asarray(t, dtype=float)))),
            'pdf': (lambda t, _r=rate: float(_r * np.exp(-_r * np.asarray(t, dtype=float)))),
            'hazard': (lambda t, _r=rate: _r),
            'mean': 1.0 / rate,
            'isExponential': True,
            'rate': rate,
        }

    from ..mam.map_analysis import map_cdf, map_pdf
    D0 = np.asarray(pair[0], dtype=float)
    D1 = np.asarray(pair[1], dtype=float)

    def ccdf(t, _D0=D0, _D1=D1):
        return float(1.0 - map_cdf(_D0, _D1, np.atleast_1d(np.asarray(t, dtype=float)))[0])

    def pdf(t, _D0=D0, _D1=D1):
        return float(map_pdf(_D0, _D1, np.atleast_1d(np.asarray(t, dtype=float)))[0])

    def hazard(t):
        # h = f/(1-F). Past the point where the ccdf underflows the hazard is
        # the asymptotic decay rate, and returning the last finite ratio is
        # better conditioned than dividing two zeros.
        c = ccdf(t)
        if c <= 1e-300:
            return rate if rate > 0 else 0.0
        return pdf(t) / c

    return {
        'ccdf': ccdf,
        'pdf': pdf,
        'hazard': hazard,
        'mean': (1.0 / rate) if rate > 0 else float('inf'),
        'isExponential': is_exp,
        'rate': rate,
    }
