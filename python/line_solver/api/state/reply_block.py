"""Synchronous-call (REPLY signal) blocked-server block in the local state.

A job of a class r with ``sn.syncreply[r] >= 0`` makes a synchronous call: it
leaves the caller station for the callee but KEEPS its server, which stays held
until the matching REPLY signal class comes back. The held servers are not
derivable from the marginal state (the job is at the callee, not here), so they
are counted per calling class in a dedicated block of the node's
local-variable vector. LDES keys the same information by job id
(Solver_ssj.pendingReplyMap); a CTMC has no job identity, so it carries counts.

Layout. The block is the LAST ``width`` columns of the node's local-variable
vector, one column per calling class that can hold a server here, in class
order. That matches MATLAB (+State/replyBlockInfo.m), where the block trails
the modulation, routing and node blocks at nvars columns 2R+1+r. The native
Python ``sn.nvars`` is an accumulated per-node width rather than the MATLAB
column-addressed layout, so the block is located from the tail of the vector;
``ctmc_ssg.space_generator_nodes`` guarantees the tail position by rotating the
reply columns behind any local variable it appends after state generation.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np


class ReplyBlockInfo(object):
    """Layout of the blocked-server block carried by one node.

    Attributes:
        classes: calling classes (0-based) that can hold a server at this node
        slot: (R,) index of each class's counter inside the local-variable
              vector, -1 when the class holds no slot here
        width: number of columns in the block
    """

    __slots__ = ('classes', 'slot', 'width')

    def __init__(self, R):
        self.classes = []
        self.slot = np.full(R, -1, dtype=int)
        self.width = 0


def reply_width(sn, ind):
    """Number of blocked-server counters node ``ind`` carries (0 when none)."""
    rb = getattr(sn, 'replyblock', None)
    if rb is None:
        return 0
    rb = np.atleast_2d(np.asarray(rb))
    if rb.size == 0 or ind >= rb.shape[0]:
        return 0
    return int(np.sum(rb[ind] > 0))


def reply_block_info(sn, ind):
    """Layout of node ``ind``'s blocked-server block (see module docstring)."""
    R = int(sn.nclasses)
    rinfo = ReplyBlockInfo(R)
    width = reply_width(sn, ind)
    if width == 0:
        return rinfo
    nvars = getattr(sn, 'nvars', None)
    if nvars is None:
        return rinfo
    nvars = np.atleast_2d(np.asarray(nvars))
    if ind >= nvars.shape[0]:
        return rinfo
    # The block sits at the tail of the local-variable vector. Before the space
    # generator has recorded the node width, nvars still reads 0; the block is
    # then the whole vector, which is what the space carries at that point.
    total = int(np.sum(nvars[ind]))
    off = max(total - width, 0)
    rb = np.atleast_2d(np.asarray(sn.replyblock))
    pos = 0
    for r in range(R):
        if rb[ind, r] > 0:
            rinfo.slot[r] = off + pos
            rinfo.classes.append(r)
            pos += 1
    rinfo.width = width
    return rinfo


def reply_blocked(sn, ind, space_var):
    """Servers held at node ``ind`` by jobs awaiting their REPLY signal.

    Returns (b, nb) with b the (rows, R) per-class counts, one row per row of
    ``space_var``, and nb the (rows,) total. Zeros when the node carries no
    block, so callers can subtract nb from the server count unconditionally.
    """
    R = int(sn.nclasses)
    sv = np.atleast_2d(np.asarray(space_var)) if space_var is not None else np.zeros((1, 0))
    nrows = max(1, sv.shape[0])
    b = np.zeros((nrows, R))
    if sv.size == 0:
        return b, np.zeros(nrows)
    rinfo = reply_block_info(sn, ind)
    if rinfo.width == 0:
        return b, np.zeros(nrows)
    for r in rinfo.classes:
        col = int(rinfo.slot[r])
        if 0 <= col < sv.shape[1]:
            b[:sv.shape[0], r] = sv[:, col]
    return b, np.sum(b, axis=1)


def reply_blocked_row(sn, ind, var_row):
    """Total servers held in a single local-variable row (scalar)."""
    _, nb = reply_blocked(sn, ind, np.atleast_2d(np.asarray(var_row)))
    return float(nb[0]) if nb.size > 0 else 0.0


def is_reply_class(sn, job_class):
    """True when ``job_class`` is a REPLY signal class."""
    st = getattr(sn, 'signaltype', None)
    if st is None or job_class >= len(st) or st[job_class] is None:
        return False
    from ...lang.classes import SignalType
    try:
        return st[job_class] == SignalType.REPLY
    except Exception:
        return False
