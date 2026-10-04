"""
Canonical (D0, D1) form for sn.proc, and a compatibility view of it.

MATLAB stores `sn.proc{i,r}` as the (D0, D1) MAP representation of the service
process (`_kb/04-networkstruct.md`). Python historically stored a per-family
DESCRIPTOR instead -- `{'rate'}` for Exp, `{'k','mu'}` for Erlang,
`{'probs','rates'}` for HyperExp, `[alpha, T]` for PH/APH -- so no reader could
treat the field uniformly and every consumer re-derived the matrices itself.

This module holds the one conversion in each direction:

- `proc_from_dist_params` builds (D0, D1) for the Markovian families, which is
  what the writer now stores.
- `proc_as_descriptor` reconstructs the historical dict/list view, so readers
  that still branch on the old shapes keep working while they are migrated.

Non-Markovian families (Gamma, Uniform, Det, Pareto, ...) keep their raw
parameter list in both codebases: `sn_nonmarkov_toph` refits the pdf from those
parameters, so there is no MAP to store until it has run.
"""

import numpy as np


def map_from_exp(rate):
    """(D0, D1) of an exponential with the given rate."""
    r = float(rate)
    return np.array([[-r]]), np.array([[r]])


def map_from_erlang(k, mu):
    """(D0, D1) of an Erlang-k whose phases each fire at rate mu."""
    k = int(k)
    mu = float(mu)
    D0 = np.zeros((k, k))
    for i in range(k):
        D0[i, i] = -mu
        if i < k - 1:
            D0[i, i + 1] = mu
    D1 = np.zeros((k, k))
    D1[k - 1, 0] = mu  # completion returns to the first phase of the next job
    return D0, D1


def map_from_hyperexp(probs, rates):
    """(D0, D1) of a hyperexponential with the given branch probs and rates."""
    p = np.asarray(probs, dtype=float).flatten()
    mu = np.asarray(rates, dtype=float).flatten()
    D0 = np.diag(-mu)
    D1 = np.outer(mu, p)  # completion re-enters branch j with probability p_j
    return D0, D1


def map_from_ph(alpha, T):
    """(D0, D1) of a phase-type law (alpha, T): D1 = (-T e) alpha."""
    T = np.asarray(T, dtype=float)
    a = np.asarray(alpha, dtype=float).flatten()
    exit_rates = -T.sum(axis=1)
    return T, np.outer(exit_rates, a)


def proc_as_descriptor(entry, procid=None):
    """The historical descriptor view of a `sn.proc` entry.

    Readers that still branch on `{'rate'}`, `{'k','mu'}`, `{'probs','rates'}`
    or `[alpha, T]` call this instead of reading the field positionally, so the
    stored form can be (D0, D1) before every reader has been migrated.

    Returns the entry unchanged when it is not a (D0, D1) pair, which covers the
    non-Markovian families and anything already in descriptor form.
    """
    if not (isinstance(entry, (list, tuple)) and len(entry) == 2):
        return entry
    D0 = np.asarray(entry[0], dtype=float)
    D1 = np.asarray(entry[1], dtype=float)
    if D0.ndim != 2 or D1.shape != D0.shape:
        return entry

    n = D0.shape[0]
    if n == 1:
        return {'rate': float(D1[0, 0])}

    off = D0 - np.diag(np.diag(D0))
    is_erlang = (np.allclose(np.diag(D0), np.diag(D0)[0])
                 and np.allclose(off, np.diag(np.full(n - 1, off[0, 1] if n > 1 else 0.0), 1))
                 and np.isclose(D1[n - 1, 0], -D0[0, 0])
                 and np.isclose(D1.sum(), -D0[0, 0]))
    if is_erlang:
        return {'k': n, 'mu': float(-D0[0, 0])}

    if np.allclose(off, 0.0):
        mu = -np.diag(D0)
        rowsum = D1.sum(axis=1)
        if np.allclose(rowsum, mu):
            p = D1[0, :] / mu[0] if mu[0] > 0 else np.full(n, 1.0 / n)
            if np.allclose(D1, np.outer(mu, p)):
                return {'probs': p, 'rates': mu}

    exit_rates = -D0.sum(axis=1)
    tot = exit_rates.sum()
    alpha = D1.sum(axis=0) / tot if tot > 0 else np.full(n, 1.0 / n)
    return [alpha, D0]


def proc_to_map(entry):
    """(D0, D1) for a `sn.proc` entry, whichever form it is stored in.

    The inverse direction of `proc_as_descriptor`: readers that want matrices
    call this and get them whether the field holds (D0, D1) already or a legacy
    descriptor. Returns (None, None) for a non-Markovian parameter list.
    """
    if isinstance(entry, dict):
        if 'k' in entry and 'mu' in entry:
            return map_from_erlang(entry['k'], entry['mu'])
        if 'rate' in entry:
            return map_from_exp(entry['rate'])
        if 'probs' in entry and 'rates' in entry:
            return map_from_hyperexp(entry['probs'], entry['rates'])
        return None, None
    # A MARKED entry is LONGER THAN TWO BLOCKS AND STILL A MAP. sn.proc holds a
    # marked arrival in the M3A layout {D0, D1, D11, ..., D1K}, so requiring
    # exactly two blocks reported "no MAP form" for every MMAP and callers read
    # the fallback of one phase -- which is how the CTMC came to pin a marked
    # Source in its first modulating phase. The (D0, D1) view of a marked entry
    # is its first two blocks; the per-mark matrices are the marking, and
    # proc_marks reads those.
    if isinstance(entry, (list, tuple)) and len(entry) >= 2:
        a = np.asarray(entry[0], dtype=float)
        b = np.asarray(entry[1], dtype=float)
        if a.ndim == 2 and b.shape == a.shape:
            return a, b          # already (D0, D1)
        if len(entry) == 2 and a.ndim <= 1 and b.ndim == 2:
            return map_from_ph(a, b)   # legacy [alpha, T]
    return None, None


def proc_to_ph(entry):
    """(alpha, T) phase-type pair for a `sn.proc` entry, in any stored form.

    The PH view of the same object `proc_to_map` returns as (D0, D1): T is D0,
    and alpha is recovered from D1 as the distribution of the phase a service
    re-enters, alpha_j = sum_i D1[i,j] / sum(D1). Readers that want a PH law
    rather than a MAP call this instead of branching on the storage shape.

    Returns (None, None) for a non-Markovian parameter list, matching what the
    per-family parsers this replaces returned when they could not identify the
    entry.
    """
    D0, D1 = proc_to_map(entry)
    if D0 is None:
        return None, None
    tot = D1.sum()
    if not tot > 0:
        n = D0.shape[0]
        return np.full(n, 1.0 / n), D0
    return D1.sum(axis=0) / tot, D0


def proc_n_phases(entry):
    """Number of phases of a `sn.proc` entry, or 1 when it has no MAP form."""
    D0, _ = proc_to_map(entry)
    return 1 if D0 is None else int(D0.shape[0])
