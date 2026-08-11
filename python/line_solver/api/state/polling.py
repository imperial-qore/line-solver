"""
Polling controller for POLLING stations.

Port of the MATLAB State.polling* helpers (matlab/src/lang/+State/polling*.m).

The polling controller is stored in the trailing columns of the local-variable
block of a station state, and holds, in order, the columns [pos, swk, ctr]; each
is materialized only when the discipline actually needs it, so that a polling
station never carries state that its dynamics cannot distinguish.

  pos  index of the buffer the server is currently at (serving) or heading to
       (switching). Materialized only when at least one switchover is
       non-immediate: while a job is in service pos always equals the class of
       that job, and while the station is empty and every switchover is
       immediate the server position is unobservable (see PARKED below).
  swk  0 when the server sits at pos, otherwise the 1-based phase of the
       switchover PH into buffer pos. Materialized only when some switchover is
       non-immediate.
  ctr  the visit budget, materialized for every discipline except EXHAUSTIVE:
         GATED        jobs of class pos admitted at the polling instant that
                      have not completed yet (the job in service counts as one
                      of them), so the visit ends when ctr reaches 0;
         KLIMITED     services still permitted in this visit, the one in
                      progress included;
         DECREMENTING the target class-pos population: the visit ends once the
                      population has dropped to ctr, i.e. one below the level
                      found at the polling instant (semi-exhaustive).

The three tangible controller configurations are therefore
  SERVING(p)   swk=0, one class-p job in the service facility;
  SWITCHING(p) swk>0, service facility empty;
  PARKED       swk=0, service facility empty, station empty. Reachable only when
               every switchover is immediate, in which case a server that
               completes a full lap without finding work would otherwise cycle
               in zero time forever. pos is then unobservable and canonical.

Immediate() switchovers are NOT represented as states: taking their rate
literally would put a ~1e8 rate in the generator (stiff, and a spurious state
per buffer). They are instead folded into the enclosing transition by
polling_next, which walks the cyclic order until a tangible state.

The block trails every other local variable (MAP server phases, round-robin
pointers, BAS marker). That ordering is load-bearing: the round-robin pointer is
located from the head of the local-variable block, so the polling columns must
come after it. Accordingly the columns here are indexed by counting back from
the end of space_var, never forward from its head.

Class indices are 0-based throughout (MATLAB uses 1-based); a service facility
holding no job is denoted by srvclass = -1.
"""

import numpy as np

from ...constants import PollingType
from ...lang.base import SchedStrategy


def _sched_val(x):
    return int(x.value) if hasattr(x, 'value') else int(x)


def _is_polling_station(sn, ind):
    """True when node ind is a station scheduled with POLLING."""
    if not sn.isstation[ind]:
        return False
    ist = int(sn.nodeToStation[ind])
    return _sched_val(sn.sched[ist]) == _sched_val(SchedStrategy.POLLING)


def _proc_disabled(sn, ist, r):
    """True when the service process of class r at station ist is absent or NaN
    (the class can never hold a job here)."""
    if sn.proc is None:
        return True
    try:
        proc_ir = sn.proc[ist][r]
    except (IndexError, KeyError, TypeError):
        return True
    if proc_ir is None:
        return True
    if isinstance(proc_ir, (list, tuple)) and len(proc_ir) > 0:
        return bool(np.any(np.isnan(np.atleast_2d(proc_ir[0]))))
    return False


def polling_info(sn, ind):
    """Derived description of the polling controller at the station of node ind,
    or None when the node is not a polling station.

    The description is a pure function of the model but is read once per state
    per synchronization while the generator is built, so it is memoized on sn.
    """
    if not _is_polling_station(sn, ind):
        return None

    cache = getattr(sn, '_pollinfo_cache', None)
    if cache is None:
        cache = {}
        sn._pollinfo_cache = cache
    if ind in cache:
        return cache[ind]

    from ..mam import map_pie

    R = int(sn.nclasses)
    ist = int(sn.nodeToStation[ind])

    pinfo = {
        'ptype': PollingType.EXHAUSTIVE,
        'pk': 1,
        'has_sw': np.zeros(R, dtype=bool),
        'ksw': np.zeros(R, dtype=int),
        'sw_d0': [None] * R,
        'sw_d1': [None] * R,
        'sw_pie': [None] * R,
    }

    # Buffers the server actually visits. A class disabled at this station can
    # never hold a job, and after_event short-circuits every event carrying it
    # (K[class] == 0), so it cannot be given a switchover to walk through: such
    # a buffer is dropped from the cyclic order rather than polled forever.
    polled = np.ones(R, dtype=bool)
    for r in range(R):
        if _proc_disabled(sn, ist, r):
            polled[r] = False
    if not np.any(polled):
        raise ValueError('Polling station %d has no class with an enabled service.' % ist)
    pinfo['polled'] = polled

    np_ind = sn.nodeparam[ind] if (sn.nodeparam is not None and ind in sn.nodeparam) else None

    def _cls_param(r):
        if not isinstance(np_ind, dict):
            return None
        v = np_ind.get(r, None)
        return v if isinstance(v, dict) else None

    for r in range(R):
        p = _cls_param(r)
        if p is None:
            continue
        pt = p.get('pollingType', None)
        if pt is not None:
            # set_polling_type writes the same discipline into every class buffer
            pinfo['ptype'] = pt
            par = p.get('pollingPar', None)
            if par:
                pinfo['pk'] = int(round(float(np.atleast_1d(par)[0])))

    # Which class buffer's switchover distribution governs the leg INTO buffer q.
    # set_switchover documents its argument as the time to switch from queue i to
    # the next one, i.e. the stored value is the cost of LEAVING buffer i, which
    # is the same leg as entering the next buffer the server visits. The
    # controller indexes its states by destination, so the lookup is shifted onto
    # the previous polled buffer here, once, rather than at every use. This is
    # also the convention of Takagi's r_i in the exact formulas of api/polling.
    polledlist = [q for q in range(R) if polled[q]]
    swof = {}
    for jj, q in enumerate(polledlist):
        swof[q] = polledlist[(jj - 1) % len(polledlist)]

    for q in polledlist:
        src = swof[q]
        p = _cls_param(src)
        if p is None:
            continue
        sw = p.get('switchoverTime', None)
        if not sw:
            continue
        # polling stores the per-buffer distribution under key 0
        dist = sw.get(0, None) if isinstance(sw, dict) else None
        if dist is None:
            continue
        if bool(dist.isImmediate()):
            continue  # zero-time switchover: folded, never a state
        repres = dist.getRepresentation()
        D0 = np.atleast_2d(np.asarray(repres[0], dtype=float))
        D1 = np.atleast_2d(np.asarray(repres[1], dtype=float))
        pinfo['has_sw'][q] = True
        pinfo['ksw'][q] = D0.shape[0]
        pinfo['sw_d0'][q] = D0
        pinfo['sw_d1'][q] = D1
        pinfo['sw_pie'][q] = np.ravel(np.asarray(map_pie(D0, D1), dtype=float))

    # Column widths. pos and swk exist only to encode SWITCHING(p); with no
    # non-immediate switchover the server is either serving (pos = class in
    # service) or parked (pos unobservable), so neither column carries
    # information. ctr exists for every discipline that bounds a visit.
    any_sw = bool(np.any(pinfo['has_sw']))
    pinfo['wpos'] = 1 if any_sw else 0
    pinfo['wswk'] = 1 if any_sw else 0
    pinfo['wctr'] = 0 if pinfo['ptype'] == PollingType.EXHAUSTIVE else 1
    pinfo['width'] = pinfo['wpos'] + pinfo['wswk'] + pinfo['wctr']

    cache[ind] = pinfo
    return pinfo


def polling_width(sn, ind):
    """Number of local-variable columns the polling controller occupies at node
    ind; 0 when the node is not a polling station."""
    pinfo = polling_info(sn, ind)
    return 0 if pinfo is None else int(pinfo['width'])


def polling_budget(pinfo, nbufq):
    """Initial value of the ctr column of a visit that starts at a buffer holding
    nbufq waiting jobs. The service facility is empty when a visit starts, so
    nbufq is the whole class population at the station."""
    ptype = pinfo['ptype']
    if ptype == PollingType.EXHAUSTIVE:
        return 0  # unused: the visit ends when the buffer drains
    if ptype == PollingType.GATED:
        return int(nbufq)  # serve exactly the jobs found at the polling instant
    if ptype == PollingType.KLIMITED:
        return int(pinfo['pk'])  # serve at most K, fewer if the buffer drains first
    if ptype == PollingType.DECREMENTING:
        return int(nbufq) - 1  # serve until the population drops one below the level found
    raise ValueError('Unsupported polling type: %s.' % str(ptype))


def polling_next(pinfo, pos, nbuf, R, arrived=False):
    """Resolve the tangible controller state a polling server reaches once it
    stops serving buffer pos. nbuf[r] is the number of class-r jobs waiting (the
    service facility is empty at this point).

    arrived selects where the cyclic walk starts. When False the server is
    LEAVING pos, so the walk starts at pos+1; when True the server has just
    ARRIVED at pos (a switchover into pos completed, or the server is parked at
    pos) and pos itself is examined first, without charging its switchover a
    second time.

    Returns (q, mode, budget):
      mode 1  start a visit to buffer q: q has work and is reached in zero time.
      mode 2  enter the switchover into buffer q: a strictly positive timer, so
              this is where the server dwells.
      mode 0  park at q: a full lap found no work and met no timed switchover,
              which can only happen when the station is empty and every
              switchover is immediate. q is canonical and unobservable.
    """
    polled = pinfo['polled']
    has_sw = pinfo['has_sw']

    if arrived and polled[pos] and nbuf[pos] > 0:
        # the switchover into pos has already been paid, so a visit starts here
        return pos, 1, polling_budget(pinfo, nbuf[pos])

    p = pos
    for _ in range(R):  # a full lap, so that the last buffer examined is pos itself
        p = (p + 1) % R
        if not polled[p]:
            continue
        if has_sw[p]:
            return p, 2, 0
        if nbuf[p] > 0:
            return p, 1, polling_budget(pinfo, nbuf[p])

    return pos, 0, 0


def polling_blocks(pinfo, srvclass, nbuf, R):
    """Every [pos, swk, ctr] controller triple compatible with a service facility
    holding a class-srvclass job (-1 when empty) and buffers holding nbuf.
    Returns an empty list when the combination is unoccupiable.

    SERVING(p): pos is pinned to srvclass; the visit budget is free within the
    bounds its discipline can have left it in. SWITCHING(q): the facility must be
    empty and q must cost time to reach; every phase of that switchover is
    reachable, and jobs may wait meanwhile, which is what makes a polling station
    non-work-conserving. PARKED: facility and station empty and no switchover
    takes time; with work waiting and only immediate switchovers the server would
    have reached it in zero time, so an idle facility with a non-empty station is
    unoccupiable and yields no rows at all. That pruning matters even when no
    column is materialized: leaving those rows in would make the generator
    reducible.
    """
    trips = []
    if srvclass >= 0:
        if not pinfo['polled'][srvclass]:
            return trips  # a buffer outside the cyclic order can hold no job
        ptype = pinfo['ptype']
        if ptype == PollingType.EXHAUSTIVE:
            ctrset = [0]
        elif ptype == PollingType.GATED:
            ctrset = range(1, int(nbuf[srvclass]) + 2)
        elif ptype == PollingType.KLIMITED:
            ctrset = range(1, int(pinfo['pk']) + 1)
        else:  # DECREMENTING
            ctrset = range(0, int(nbuf[srvclass]) + 1)
        for ctr in ctrset:
            trips.append((srvclass, 0, int(ctr)))
    else:
        for q in np.where(pinfo['has_sw'])[0]:
            for swk in range(1, int(pinfo['ksw'][q]) + 1):
                trips.append((int(q), swk, 0))
        if not np.any(pinfo['has_sw']) and int(np.sum(nbuf)) == 0:
            trips.append((int(np.where(pinfo['polled'])[0][0]), 0, 0))  # parked, pos canonical
    return trips


def polling_project(pinfo, trips):
    """Project full [pos, swk, ctr] triples onto the columns polling_info
    materializes. The elided columns are reconstructible from the rest of the
    state (see polling_get), so keeping them would split each state into copies
    no observation can tell apart."""
    keep = np.array([pinfo['wpos'], pinfo['wswk'], pinfo['wctr']], dtype=bool)
    arr = np.atleast_2d(np.asarray(trips, dtype=float))
    if arr.size == 0:
        return np.zeros((0, int(pinfo['width'])))
    return arr[:, keep]


def _base(pinfo, nvar_cols):
    """First column of the polling block within a space_var row of nvar_cols
    columns. The block trails every other local variable, so it is located by
    counting back from the end."""
    return int(nvar_cols) - int(pinfo['width'])


def polling_get(pinfo, space_var_row, srvclass):
    """Read the polling controller out of the local-variable columns of a single
    state row. srvclass is the class in the service facility, -1 when empty.

    Columns that polling_info elides are reconstructed here, so callers always
    see a complete controller: pos falls back to the class in service (or the
    canonical first polled buffer when parked), swk to 0 (no switchover takes
    time), ctr to 0 (EXHAUSTIVE bounds a visit by the buffer draining).
    """
    row = np.ravel(np.asarray(space_var_row, dtype=float))
    b = _base(pinfo, row.size)
    c = b
    if pinfo['wpos']:
        pos = int(round(row[c]))
        c += 1
    elif srvclass >= 0:
        pos = int(srvclass)
    else:
        pos = int(np.where(pinfo['polled'])[0][0])
    if pinfo['wswk']:
        swk = int(round(row[c]))
        c += 1
    else:
        swk = 0
    if pinfo['wctr']:
        ctr = int(round(row[c]))
    else:
        ctr = 0
    return pos, swk, ctr


def polling_set(pinfo, space_var_row, pos, swk, ctr):
    """Write the polling controller into the local-variable columns of a state
    row, returning a new row. Columns polling_info elides are dropped."""
    row = np.ravel(np.asarray(space_var_row, dtype=float)).copy()
    b = _base(pinfo, row.size)
    c = b
    if pinfo['wpos']:
        row[c] = pos
        c += 1
    if pinfo['wswk']:
        row[c] = swk
        c += 1
    if pinfo['wctr']:
        row[c] = ctr
    return row


def polling_land(pinfo, q, mode, budget, space_buf, space_srv, space_var,
                 K, Ks, pie, ist, R):
    """Materialize the state rows a polling server lands in after polling_next
    resolved (q, mode, budget), with the probability of each. The three space_*
    inputs are single rows describing the station at the instant the decision is
    taken, i.e. with the completed job (if any) already removed from the server.

    probs splits a landing across the entry phases of a phase-type: which phase a
    service or a switchover starts in is a random choice, so one decision yields
    one row per entry phase. Callers fold probs into the transition rate rather
    than into outprob, since the branching happens when the active event fires.
    """
    from .after_event_station import _get_pie

    buf = np.ravel(np.asarray(space_buf, dtype=float))
    srv = np.ravel(np.asarray(space_srv, dtype=float))
    var = np.ravel(np.asarray(space_var, dtype=float))
    rows = []
    probs = []

    if mode == 1:
        # start or continue a visit at q: pull a waiting class-q job into service
        pentry = _get_pie(pie, ist, q, K)
        for kentry in range(int(K[q])):
            if pentry[kentry] <= 0:
                continue
            buf_k = buf.copy()
            buf_k[q] -= 1
            srv_k = srv.copy()
            srv_k[int(Ks[q]) + kentry] += 1
            var_k = polling_set(pinfo, var, q, 0, budget)
            rows.append(np.concatenate([buf_k, srv_k, var_k]))
            probs.append(float(pentry[kentry]))
    elif mode == 2:
        # enter the switchover into q: the facility stays empty while walking
        swpie = pinfo['sw_pie'][q]
        for kentry in range(int(pinfo['ksw'][q])):
            if swpie[kentry] <= 0:
                continue
            var_k = polling_set(pinfo, var, q, kentry + 1, 0)
            rows.append(np.concatenate([buf, srv, var_k]))
            probs.append(float(swpie[kentry]))
    else:
        # park: held until the next arrival, see polling_next
        var_k = polling_set(pinfo, var, q, 0, 0)
        rows.append(np.concatenate([buf, srv, var_k]))
        probs.append(1.0)

    return rows, probs


def _srvclass_of(srv_row, K, Ks, R):
    """Class occupying the single service facility of a polling station, -1 when
    it is empty."""
    srv = np.ravel(np.asarray(srv_row, dtype=float))
    for r in range(R):
        if float(np.sum(srv[int(Ks[r]):int(Ks[r]) + int(K[r])])) > 0:
            return r
    return -1


def polling_space(sn, ind, space, K, Ks):
    """Append the polling controller columns to the rows of space, which must
    hold the [buffer, server, other-local-variable] layout of a polling station.
    Rows are expanded into one row per controller configuration the discipline
    can occupy, and rows for which no configuration exists are dropped.

    Enumerating the controller per row rather than as a blind cartesian product
    is what keeps the state space tight and the chain irreducible: pos is pinned
    to the class in service, a switchover excludes a busy service facility, and a
    park excludes a non-empty station. A cartesian product would instead admit
    states such as "serving buffer 1 while a class-2 job occupies the server",
    which no transition can reach or leave consistently with the marginals.
    """
    pinfo = polling_info(sn, ind)
    space = np.atleast_2d(np.asarray(space, dtype=float))
    if pinfo is None or space.size == 0:
        return space

    R = int(sn.nclasses)
    sumK = int(np.sum(K))
    out = []
    for row in range(space.shape[0]):
        nbuf = space[row, :R]
        srv = space[row, R:R + sumK]
        srvclass = _srvclass_of(srv, K, Ks, R)
        blocks = polling_project(pinfo, polling_blocks(pinfo, srvclass, nbuf, R))
        for b in range(blocks.shape[0]):
            out.append(np.concatenate([space[row, :], blocks[b, :]]))
    if not out:
        return np.zeros((0, space.shape[1] + int(pinfo['width'])))
    return np.array(out, dtype=float)


def polling_init(sn, ind, nbuf, srvclass):
    """Canonical initial value of the polling controller columns for a station
    whose buffers hold nbuf and whose facility holds a class-srvclass job (-1
    when empty). Returns a 1-D array of pinfo['width'] entries (possibly empty).

    The row is always one polling_blocks also enumerates, so the initial state is
    a member of the generated state space.
    """
    pinfo = polling_info(sn, ind)
    if pinfo is None:
        return np.zeros(0)
    R = int(sn.nclasses)
    nbuf = np.ravel(np.asarray(nbuf, dtype=float))

    if srvclass >= 0:
        # A visit to srvclass is under way. Take it to have started with the
        # whole class present, which is the state reached by a server that has
        # just arrived at a buffer holding nbuf[srvclass]+1 jobs and pulled one
        # of them into service.
        pos = int(srvclass)
        swk = 0
        ctr = polling_budget(pinfo, nbuf[srvclass] + 1)
    else:
        # The facility is empty: let the server walk from a canonical position
        # and settle wherever the discipline puts it. With work waiting this
        # cannot return a visit, because initialization always fills the single
        # server of a non-empty station, so only a switchover or a park is
        # reachable here.
        start = int(np.where(pinfo['polled'])[0][0])
        q, mode, _ = polling_next(pinfo, start, nbuf, R, arrived=True)
        pos = q
        ctr = 0
        if mode == 2:
            swk = int(np.where(pinfo['sw_pie'][q] > 0)[0][0]) + 1
        else:
            swk = 0

    return np.ravel(polling_project(pinfo, [(pos, swk, ctr)]))
