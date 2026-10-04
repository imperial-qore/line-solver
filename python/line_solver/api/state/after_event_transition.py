"""
Transition node event handler for afterEvent dispatch.

Handles PHASE, ENABLE, and FIRE events at Transition nodes (SPN).
In the CTMC context, only PHASE events produce local state changes;
ENABLE/FIRE are global events handled elsewhere.

Port from JAR AfterEventTransition.java.
"""

import numpy as np
from ...constants import EventType


def after_event_transition(sn, ind, event, job_class, inspace,
                           K, Ks, space_buf, space_srv, space_var,
                           space_fired=None):
    """
    Handle events at a Transition node.

    Args:
        sn: NetworkStruct
        ind: Node index (0-based)
        event: EventType (PHASE, ENABLE, or FIRE)
        job_class: Mode index (0-based) for PHASE events
        inspace: Full input state, shape (n_rows, n_cols)
        K: Phase counts per mode, shape (nmodes,)
        Ks: Phase shift per mode, shape (nmodes,)
        space_buf: Buffer state (idle server counts per mode)
        space_srv: Server state (enabled phase counts)
        space_var: Variable state

    Returns:
        Tuple of (outspace, outrate, outprob)
    """
    inspace = np.atleast_2d(inspace)

    if event in (EventType.ENABLE, EventType.FIRE):
        # Global events: no local state change
        outspace = inspace.copy()
        outrate = np.zeros((outspace.shape[0], 1))
        outprob = np.ones((outspace.shape[0], 1))
        return outspace, outrate, outprob

    if event == EventType.PHASE:
        mode = job_class  # In a Transition, jobClass is interpreted as the mode

        # Bounds check
        if mode < 0 or mode >= len(K):
            return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

        # Single-phase modes have no phase transitions
        if K[mode] <= 1:
            return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))

        # Compute marginals via toMarginal
        from .marginal import toMarginal
        ni, nir, sir, kir = toMarginal(sn, ind, inspace, K, Ks,
                                        space_buf=space_buf,
                                        space_srv=space_srv,
                                        space_var=space_var)
        nir = np.atleast_2d(nir)
        kir = np.atleast_3d(kir) if kir.ndim < 3 else kir

        if nir[0, mode] <= 0:
            # No enabled servers in this mode
            outspace = inspace.copy()
            outrate = np.zeros((outspace.shape[0], 1))
            outprob = np.ones((outspace.shape[0], 1))
            return outspace, outrate, outprob

        # Get firing process D0 matrix for this mode
        nparam = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
        firing_proc_d0 = None
        if nparam is not None:
            fp = nparam.get('firingproc', None) if isinstance(nparam, dict) else getattr(nparam, 'firingproc', None)
            if fp is not None and mode < len(fp):
                proc_cell = fp[mode]
                if proc_cell is not None and len(proc_cell) > 0:
                    firing_proc_d0 = np.atleast_2d(proc_cell[0])

        out_states = []
        out_rates = []
        out_probs = []

        n_rows = space_srv.shape[0]

        for k in range(int(K[mode])):
            col_k = int(Ks[mode]) + k
            # Find enabled rows: those with spaceSrv[row, col_k] > 0
            en = space_srv[:, col_k] > 0
            if not np.any(en):
                continue

            enabled_rows = np.where(en)[0]

            for kdest in range(int(K[mode])):
                if kdest == k:
                    continue

                col_kdest = int(Ks[mode]) + kdest

                # Extract enabled rows
                space_srv_k = space_srv[enabled_rows].copy()
                space_buf_k = space_buf[enabled_rows].copy() if space_buf.size > 0 else np.zeros((len(enabled_rows), 0))
                space_var_k = space_var[enabled_rows].copy() if space_var.size > 0 else np.zeros((len(enabled_rows), 0))
                if space_fired is not None and space_fired.size > 0:
                    space_fired_k = np.atleast_2d(space_fired)[enabled_rows].copy()
                else:
                    space_fired_k = np.zeros((len(enabled_rows), 0))

                # Move job from phase k to phase kdest
                space_srv_k[:, col_k] -= 1
                space_srv_k[:, col_kdest] += 1

                # Compute rates per row. MATLAB afterEventTransition.m line 38
                # uses rate = D0(k,kdest) * kir(:,mode,k) and outrate = nir(mode) .* rate.
                # However, this double-multiplies by per-mode count: nir(mode)
                # already equals sum over k of kir(:,mode,k). For multi-server
                # operation the proper rate is D0(k,kdest)*kir(:,mode,k) only.
                rates = np.zeros(len(enabled_rows))
                for i, orig_row in enumerate(enabled_rows):
                    kir_val = kir[orig_row, mode, k] if kir.ndim == 3 else kir[orig_row, k]

                    if firing_proc_d0 is not None:
                        fp_rate = firing_proc_d0[k, kdest]
                    else:
                        fp_rate = 1.0

                    rates[i] = fp_rate * kir_val

                # Build output state rows in [buf | srv | fired | var] layout
                # to match sn.space[isf_transition].
                parts = [space_buf_k, space_srv_k, space_fired_k, space_var_k]
                state_rows = np.hstack([p for p in parts if p.size > 0])

                for i in range(len(enabled_rows)):
                    out_states.append(state_rows[i])
                    out_rates.append(rates[i])
                    out_probs.append(1.0)

        if len(out_states) == 0:
            # No valid transitions found
            outspace = inspace.copy()
            outrate = np.zeros((outspace.shape[0], 1))
            outprob = np.ones((outspace.shape[0], 1))
            return outspace, outrate, outprob

        outspace = np.array(out_states)
        outrate = np.array(out_rates).reshape(-1, 1)
        outprob = np.array(out_probs).reshape(-1, 1)
        return outspace, outrate, outprob

    # Unknown event
    return np.zeros((0, 0)), np.zeros((0, 0)), np.zeros((0, 0))
