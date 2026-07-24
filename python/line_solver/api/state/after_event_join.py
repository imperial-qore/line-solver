"""
afterEvent handler for Join stations on FJ-augmented structs.

Port of matlab/src/lang/+State/afterEventJoin.m.
"""

import numpy as np
from ...constants import EventType, GlobalConstants


def after_event_join(sn, ind, inspace, event, job_class, is_simulation=False):
    """
    Handle afterEvent logic for Join stations on FJ-augmented structs
    (see ModelAdapter.fjtag). The join state is a plain per-class count
    vector [n_1..n_R] of buffered jobs: auxiliary-class entries count the
    sibling tasks waiting for synchronization.

    ARV (passive): buffer the arriving job/sibling.
    DEP in an original class r (active, immediate): enabled when a plain
    class-r job is buffered, or when some tag t has the full required
    sibling multiset present (all B branches, identity matching by tag).
    The firing consumes the siblings of the LOWEST complete tag and lets
    one class-r job depart.

    Returns (outspace, outrate, outprob).
    """
    R = int(sn.nclasses)
    inspace = np.atleast_2d(np.array(inspace, dtype=float))
    n = inspace[:, -R:].copy()  # per-class counts (rows x R)

    outspace = None
    outrate = None
    outprob = np.ones((1, 1))

    if event == EventType.ARV:
        n[:, job_class] = n[:, job_class] + 1
        outspace = n
        outrate = -1.0 * np.ones((outspace.shape[0], 1))  # passive, rate unspecified
        outprob = np.ones((outspace.shape[0], 1))
    elif event == EventType.DEP:
        fjp = None
        if sn.nodeparam is not None and ind in sn.nodeparam and isinstance(sn.nodeparam[ind], dict):
            fjp = sn.nodeparam[ind].get('fj', None)

        # auxiliary sibling classes never depart individually: they are
        # consumed only by the original-class join firing below
        isaux = False
        if fjp is not None:
            for rr in np.asarray(fjp['origclasses']).ravel():
                auxm = fjp['auxmatrix'][int(rr)]
                if auxm is not None and np.any(np.asarray(auxm) == job_class):
                    isaux = True
                    break
        if isaux:
            return None, None, np.ones((1, 1))

        rows_out = []
        for row in range(n.shape[0]):
            nrow = n[row, :].copy()
            if nrow[job_class] > 0:
                # plain (non-forked) job buffered in the departing class
                nrow[job_class] = nrow[job_class] - 1
                rows_out.append(nrow)
            elif fjp is not None and job_class < len(fjp['auxmatrix']) and fjp['auxmatrix'][job_class] is not None:
                auxm = np.asarray(fjp['auxmatrix'][job_class])  # B x T
                req = np.asarray(fjp['required'][job_class]).ravel()  # B
                T = auxm.shape[1]
                tstar = -1
                for t in range(T):
                    if np.all(nrow[auxm[:, t].astype(int)] >= req):
                        tstar = t  # lowest complete tag, canonical consumption
                        break
                if tstar >= 0:
                    nrow[auxm[:, tstar].astype(int)] = nrow[auxm[:, tstar].astype(int)] - req
                    rows_out.append(nrow)

        if rows_out:
            outspace = np.array(rows_out, dtype=float)
            outrate = GlobalConstants.Immediate * np.ones((outspace.shape[0], 1))
            outprob = np.ones((outspace.shape[0], 1))
        else:
            outspace = None
            outrate = None
            outprob = np.ones((1, 1))

    if is_simulation and outspace is not None and outspace.shape[0] > 1:
        tot_rate = np.sum(outrate)
        cum_rate = np.cumsum(outrate) / tot_rate
        firing_ctr = int(np.sum(np.random.rand() > cum_rate.ravel()))
        firing_ctr = min(firing_ctr, outspace.shape[0] - 1)
        outspace = outspace[firing_ctr:firing_ctr + 1, :]
        outrate = np.array([[tot_rate]])
        outprob = outprob[firing_ctr:firing_ctr + 1, :]

    return outspace, outrate, outprob
