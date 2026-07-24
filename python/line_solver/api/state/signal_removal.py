"""G-network signal arrival at a station (port of MATLAB +State/afterEventStationSignal.m).

A signal class never joins a station: it acts on the jobs already there and is
annihilated. CATASTROPHE empties the station, ignoring the removal count
distribution. NEGATIVE removes a batch of jobs, where

  - victims are the signal's target class when one is declared
    (sn.signaltarget >= 0, i.e. forJobClass was used) and otherwise every
    non-signal class -- the untargeted case is the classic Gelenbe negative
    customer and agrees with SolverMAM and SolverLDES, which are class-agnostic;
  - the batch size follows sn.signalremdist, clipped at the eligible population
    so that the pmf tail P(B >= n) lumps onto "remove all n";
  - sn.signalrempolicy picks the victim: FCFS the oldest waiting job, LCFS the
    newest, RANDOM uniformly over waiting and in-service jobs alike. FCFS/LCFS
    only reach into the servers once the waiting line is empty (two-tier, as in
    SolverLDES).

The station state records no arrival-time order for in-service jobs, so when a
policy has to reach into the servers the victim is drawn uniformly across the
occupied phases (exact whenever at most one job of the class is in service).
"""

import numpy as np

from ...lang.base import SchedStrategy

# Buffer layouts. Only the ordered layouts record arrival order.
_ORDERED_SCHED = (SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.LCFS)
_COUNT_SCHED = (SchedStrategy.SIRO, SchedStrategy.SEPT, SchedStrategy.LEPT)


def _sched_attr(name):
    return getattr(SchedStrategy, name, None)


def _pair_sched():
    names = ('FCFSPR', 'FCFSPI', 'FCFSPRPRIO', 'FCFSPIPRIO',
             'LCFSPR', 'LCFSPI', 'LCFSPRPRIO', 'LCFSPIPRIO')
    return tuple(s for s in (_sched_attr(n) for n in names) if s is not None)


def _ordered_sched():
    extra = tuple(s for s in (_sched_attr('LCFSPRIO'),) if s is not None)
    return _ORDERED_SCHED + extra


def is_catastrophe_signal(sn, job_class):
    """True if the signal empties a station on arrival.

    signaltype is consulted alongside the iscatastrophe flag so the two
    encodings of a catastrophe cannot disagree (SolverMAM applies the same test).
    """
    isc = getattr(sn, 'iscatastrophe', None)
    if isc is not None and job_class < len(isc) and bool(isc[job_class]):
        return True
    st = getattr(sn, 'signaltype', None)
    if st is not None and job_class < len(st) and st[job_class] is not None:
        from ...lang.classes import SignalType
        try:
            return st[job_class] == SignalType.CATASTROPHE
        except Exception:
            return False
    return False


def signal_batch_pmf(sn, job_class, ntot):
    """Batch-size distribution, clipped at the ntot eligible jobs.

    Without a removal distribution a signal removes exactly one job. With one,
    a batch larger than the eligible population empties it rather than driving
    the queue negative, so P(B >= ntot) lumps onto "remove ntot". This matches
    the min(B, n) clipping in SolverLDES and the tail term SolverMAM puts on the
    empty state.
    """
    remdist = getattr(sn, 'signalremdist', None)
    if remdist is None or job_class >= len(remdist) or remdist[job_class] is None:
        return [1], [1.0]
    dist = remdist[job_class]

    head = []
    for b in range(int(ntot)):
        head.append(float(dist.evalPMF(b)))
    tail = max(0.0, 1.0 - float(np.sum(head)))
    kvals = list(range(int(ntot))) + [int(ntot)]
    kprobs = head + [tail]

    keep = [(k, p) for k, p in zip(kvals, kprobs) if p > 0]
    if not keep:
        return [1], [1.0]
    kvals = [k for k, _ in keep]
    kprobs = [p for _, p in keep]
    total = float(np.sum(kprobs))
    kprobs = [p / total for p in kprobs]
    return kvals, kprobs


def _merge_states(space, prob):
    if len(space) == 0:
        return space, prob
    arr = np.asarray(space, dtype=float)
    uniq, inv = np.unique(arr, axis=0, return_inverse=True)
    p = np.zeros(uniq.shape[0])
    for i, idx in enumerate(np.ravel(inv)):
        p[idx] += prob[i]
    return [uniq[i, :] for i in range(uniq.shape[0])], list(p)


def _waiting_victims(sched, buf, tgtclasses):
    """Eligible waiting jobs as (positions, classes, multiplicities, ordered, pair)."""
    pos, cls, weight = [], [], []
    is_ordered = False
    is_pair = False
    if buf.size == 0:
        return pos, cls, weight, is_ordered, is_pair
    if sched in _ordered_sched():
        is_ordered = True
        for c in range(len(buf)):
            v = int(round(buf[c]))
            if v > 0 and (v - 1) in tgtclasses:
                pos.append(c)
                cls.append(v - 1)
                weight.append(1.0)
    elif sched in _pair_sched():
        is_ordered = True
        is_pair = True
        for c in range(0, len(buf) - 1, 2):
            v = int(round(buf[c]))
            if v > 0 and (v - 1) in tgtclasses:
                pos.append(c)
                cls.append(v - 1)
                weight.append(1.0)
    elif sched in _COUNT_SCHED:
        # per-class counts: each eligible class is one victim kind carrying as
        # many interchangeable jobs as its count
        for r in tgtclasses:
            if r < len(buf) and buf[r] > 0:
                pos.append(r)
                cls.append(r)
                weight.append(float(buf[r]))
    # PS/INF/DPS/GPS/LPS and the source hold no waiting jobs
    return pos, cls, weight, is_ordered, is_pair


def _in_service_victims(srv, tgtclasses, K, Ks):
    cls, phase, count = [], [], []
    for r in tgtclasses:
        for p in range(int(K[r])):
            n = srv[int(Ks[r]) + p]
            if n > 0:
                cls.append(r)
                phase.append(p)
                count.append(float(n))
    return cls, phase, count


def _drop_waiting(buf, pos, is_ordered, is_pair, cls):
    """Remove a waiting job, keeping the layout invariant: empty slots pad the
    left, the head of line stays rightmost."""
    buf = np.array(buf, dtype=float)
    if is_pair:
        buf = np.delete(buf, [pos, pos + 1])
        buf = np.concatenate([[0.0, 0.0], buf])
    elif is_ordered:
        buf = np.delete(buf, pos)
        buf = np.concatenate([[0.0], buf])
    else:
        buf[cls] -= 1  # per-class count buffer
    return buf


def _drop_in_service(sched, buf, srv, cls, phase, Ks, S_ist):
    """Remove an in-service job and, at a station that keeps a waiting line,
    pull the head of line into the freed server."""
    buf = np.array(buf, dtype=float)
    srv = np.array(srv, dtype=float)
    srv[int(Ks[cls]) + phase] -= 1
    if buf.size == 0 or np.sum(srv) >= S_ist:
        return buf, srv
    if sched in _ordered_sched():
        headpos = next((c for c in range(len(buf) - 1, -1, -1) if buf[c] > 0), -1)
        if headpos >= 0:
            promo = int(round(buf[headpos]))
            # the slot is vacated, not blanked in place: the waiting line stays
            # right-aligned with the empty slots padding the left
            buf = np.delete(buf, headpos)
            buf = np.concatenate([[0.0], buf])
            srv[int(Ks[promo - 1])] += 1
    elif sched in _pair_sched():
        headpos = -1
        for c in range(0, len(buf) - 1, 2):
            if buf[c] > 0:
                headpos = c
        if headpos >= 0:
            promo = int(round(buf[headpos]))
            promophase = int(round(buf[headpos + 1]))
            if promophase < 1:
                promophase = 1
            buf = np.delete(buf, [headpos, headpos + 1])
            buf = np.concatenate([[0.0, 0.0], buf])
            srv[int(Ks[promo - 1]) + promophase - 1] += 1  # resumes at its stored phase
    elif sched in _COUNT_SCHED:
        promo = next((c for c in range(len(buf)) if buf[c] > 0), -1)
        if promo >= 0:
            buf[promo] -= 1
            srv[int(Ks[promo])] += 1
    return buf, srv


def _remove_one(sched, buf, srv, var, tgtclasses, policy, K, Ks, S_ist):
    from ...lang.classes import RemovalPolicy

    outspace, outprob = [], []
    pos, cls, weight, is_ordered, is_pair = _waiting_victims(sched, buf, tgtclasses)
    scls, sphase, scount = _in_service_victims(srv, tgtclasses, K, Ks)
    nwait = float(np.sum(weight)) if weight else 0.0
    nsrv = float(np.sum(scount)) if scount else 0.0
    if nwait == 0 and nsrv == 0:
        return outspace, outprob

    # FCFS/LCFS rank the waiting line by age, which only an ordered buffer
    # records; a per-class count buffer carries no age, so an age-based policy
    # degenerates to a uniform draw over the waiting jobs.
    age_ordered = is_ordered and policy in (RemovalPolicy.FCFS, RemovalPolicy.LCFS)
    if age_ordered and nwait > 0:
        if policy == RemovalPolicy.FCFS:
            pick = int(np.argmax(pos))  # head of line: the last occupied slot
        else:
            pick = int(np.argmin(pos))  # most recent arrival: the first occupied slot
        b2 = _drop_waiting(buf, pos[pick], is_ordered, is_pair, cls[pick])
        outspace.append(np.concatenate([b2, srv, var]))
        outprob.append(1.0)
        return outspace, outprob

    if policy == RemovalPolicy.RANDOM:
        total = nwait + nsrv  # uniform over waiting and in-service jobs alike
    else:
        total = nwait if nwait > 0 else nsrv  # drain the waiting line first

    if nwait > 0:
        for w in range(len(pos)):
            b2 = _drop_waiting(buf, pos[w], is_ordered, is_pair, cls[w])
            outspace.append(np.concatenate([b2, srv, var]))
            outprob.append(weight[w] / total)
    if policy == RemovalPolicy.RANDOM or nwait == 0:
        for j in range(len(scls)):
            b2, s2 = _drop_in_service(sched, buf, srv, scls[j], sphase[j], Ks, S_ist)
            outspace.append(np.concatenate([b2, s2, var]))
            outprob.append(scount[j] / total)
    return _merge_states(outspace, outprob)


def _remove_batch(sched, buf, srv, var, k, tgtclasses, policy, K, Ks, S_ist):
    """Remove k jobs one at a time; sequential uniform draws without replacement
    reproduce a uniform choice of the removed subset."""
    nb = len(buf)
    ns = len(srv)
    space = [np.concatenate([buf, srv, var])]
    prob = [1.0]
    for _ in range(int(k)):
        nextspace, nextprob = [], []
        for row in range(len(space)):
            st = space[row]
            b = st[:nb]
            s = st[nb:nb + ns]
            v = st[nb + ns:]
            sp, pr = _remove_one(sched, b, s, v, tgtclasses, policy, K, Ks, S_ist)
            if not sp:
                # nothing left to remove: the state is already drained
                nextspace.append(st)
                nextprob.append(prob[row])
            else:
                for j in range(len(sp)):
                    nextspace.append(sp[j])
                    nextprob.append(prob[row] * pr[j])
        space, prob = _merge_states(nextspace, nextprob)
    return space, prob


def handle_signal_arrival(sn, ind, ist, inspace, job_class, sched, K, Ks, S,
                          space_buf, space_srv, space_var, is_simulation=False):
    """Passive arrival of a G-network signal class at station ist.

    Which job a signal removes is in general a random choice, so the generator
    needs every destination state and its probability. A simulation instead walks
    one sample path, so under is_simulation the set is collapsed to a single
    successor drawn from outprob (the rate stays -1, as for the other passive
    actions).
    """
    from .marginal import toMarginal
    from ...lang.classes import RemovalPolicy

    _, s_nir, s_sir, _ = toMarginal(sn, ind, inspace)
    buf = np.ravel(space_buf).astype(float)
    srv = np.ravel(space_srv).astype(float)
    var = np.ravel(space_var).astype(float) if space_var.size > 0 else np.array([])

    if is_catastrophe_signal(sn, job_class):
        out = np.concatenate([np.zeros(len(buf)), np.zeros(len(srv)), var])
        return out.reshape(1, -1), np.array([[-1.0]]), np.array([[1.0]])

    tgt = -1
    st = getattr(sn, 'signaltarget', None)
    if st is not None and job_class < len(st):
        tgt = int(st[job_class])
    if tgt >= 0:
        tgtclasses = [tgt]
    else:
        tgtclasses = [r for r in range(len(sn.issignal)) if not bool(sn.issignal[r])]
    tgtclasses = [r for r in tgtclasses if float(s_nir[0, r]) > 0]
    ntot = int(round(sum(float(s_nir[0, r]) for r in tgtclasses)))
    if not tgtclasses or ntot <= 0:
        out = np.concatenate([buf, srv, var])  # no victim: the signal vanishes
        return out.reshape(1, -1), np.array([[-1.0]]), np.array([[1.0]])

    kvals, kprobs = signal_batch_pmf(sn, job_class, ntot)

    policy = RemovalPolicy.RANDOM
    pol = getattr(sn, 'signalrempolicy', None)
    if pol is not None and job_class < len(pol) and pol[job_class] is not None:
        policy = pol[job_class]

    S_ist = float(S[ist]) if np.ndim(S) > 0 else float(S)
    outspace, outprob = [], []
    for k, pk in zip(kvals, kprobs):
        if pk <= 0:
            continue
        if k <= 0:
            outspace.append(np.concatenate([buf, srv, var]))
            outprob.append(pk)
            continue
        sp, pr = _remove_batch(sched, buf, srv, var, k, tgtclasses, policy, K, Ks, S_ist)
        for j in range(len(sp)):
            outspace.append(sp[j])
            outprob.append(pk * pr[j])

    outspace, outprob = _merge_states(outspace, outprob)
    if is_simulation and len(outspace) > 1:
        cum = np.cumsum(outprob) / float(np.sum(outprob))
        rnd = np.random.rand()
        fc = 1 + max([-1] + [i for i in range(len(cum)) if rnd > cum[i]])
        return np.ravel(outspace[fc]).reshape(1, -1), np.array([[-1.0]]), np.array([[1.0]])
    out = np.vstack([np.ravel(s) for s in outspace])
    orate = -1.0 * np.ones((out.shape[0], 1))
    oprob = np.asarray(outprob, dtype=float).reshape(-1, 1)
    return out, orate, oprob
