"""Multiserver handling policy for SolverNC.

Twin of ``matlab/src/solvers/NC/nc_multiserver_policy.m`` and
``nc_lld_from_nservers.m``. See ``_kb/06-solver-catalog.md`` (NC section,
multiserver handling).
"""

import numpy as np

__all__ = ['nc_multiserver_policy', 'nc_lld_from_nservers']


def nc_multiserver_policy(options, warn=None):
    """Resolve ``options.config.multiserver`` into SolverNC's multiserver handling.

    SolverNC represents a finite multiserver station either by Seidmann's
    approximation (demand ``L/c`` plus a delay ``L(c-1)/c``) or by the exact
    load-dependent lattice ``mu(n)=min(n,c)``, which routes the model to the
    load-dependent solver. Returns one of:

    ``'default'``
        the historical dispatch: Seidmann on method ``default``, the lattice on
        ``exact``/``is``/``panald``, and the lattice on the 2-station
        Delay+multiserver topology.
    ``'seidmann'``
        Seidmann everywhere, including on ``exact``.
    ``'lld'``
        the lattice everywhere it is admissible, including on ``default``.

    ``config.multiserver`` belongs to the general solver options, shared with
    SolverMVA, which implements approximations SolverNC has no counterpart for
    (``softmin``, ``conway``, ``krzesinski``, ``suri``, ``erlang``). Those warn
    and fall back to ``'default'`` rather than erroring, because one options
    object is commonly reused across solvers.
    """
    requested = None
    config = None
    # options is a mapping in the native python solvers and an attribute-style
    # object in the wrappers; accept either
    if hasattr(options, 'get'):
        try:
            config = options.get('config', None)
        except TypeError:
            config = None
    if config is None:
        config = getattr(options, 'config', None)
    if config is not None:
        if hasattr(config, 'get'):
            requested = config.get('multiserver', None)
        else:
            requested = getattr(config, 'multiserver', None)
    if requested is None or requested == '':
        return 'default'

    requested = str(requested).lower()
    if requested == 'default':
        return 'default'
    if requested == 'seidmann':
        return 'seidmann'
    if requested in ('lld', 'exact', 'loaddep', 'load-dependent'):
        return 'lld'
    if warn is not None:
        warn("SolverNC does not implement config.multiserver='%s' (it is a SolverMVA "
             "approximation); using 'default'. SolverNC accepts 'default', 'seidmann' "
             "and 'lld'." % requested)
    return 'default'


def nc_lld_from_nservers(sn, nservers, lattice_max=None):
    """Exact ``mu(n)=min(n,c)`` lattice for a closed model's multiserver stations.

    Returns ``None`` when the conversion does not apply -- no finite multiserver
    station, an open or mixed model, an lldscaling already installed, or a
    per-chain population lattice above ``lattice_max`` -- in which case the
    caller keeps Seidmann's approximation.
    """
    existing = getattr(sn, 'lldscaling', None)
    if existing is not None and getattr(existing, 'size', 0) > 0:
        return None
    if not any(s > 1 and np.isfinite(s) for s in nservers):
        return None

    njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(sn.nclasses)
    if any(np.isinf(njobs)):
        return None
    Nt = int(np.sum(njobs[np.isfinite(njobs)]))
    if Nt < 1:
        return None

    if lattice_max is not None:
        chains = np.asarray(sn.chains, dtype=float) if sn.chains is not None else np.array([])
        lattice = 1.0
        if chains.ndim == 2 and chains.shape[0] > 0:
            for c in range(chains.shape[0]):
                popc = float(np.sum([njobs[r] for r in np.where(chains[c, :] > 0)[0]
                                     if r < njobs.size and np.isfinite(njobs[r])]))
                lattice *= (1.0 + popc)
        else:
            for v in njobs:
                if np.isfinite(v):
                    lattice *= (1.0 + float(v))
        if lattice > lattice_max:
            return None

    lldscaling = np.ones((sn.nstations, Nt))
    for i in range(sn.nstations):
        if nservers[i] > 1 and np.isfinite(nservers[i]):
            for j in range(Nt):
                lldscaling[i, j] = min(j + 1, nservers[i])
    return lldscaling
