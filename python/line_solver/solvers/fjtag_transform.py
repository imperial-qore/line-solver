"""Fork-join TAG AUGMENTATION, the transform/lift pair shared by CTMC and SSA.

``ModelAdapter.fjtag`` is the OTHER fork-join route: where ``mmt``/``ht`` drive
an outer fixed point for MVA, NC and FLD (see :mod:`fork_join_driver`), the tag
augmentation is EXACT and single pass. It rewrites the model so each sibling
branch carries its own auxiliary class, the solver runs unchanged on that
struct, and the auxiliary columns are folded back at the end.

This is deliberately NOT built on the ``fork_join_driver`` shape. That driver
owns a loop and reaches the inner solve through a seam because MMT re-solves a
transformed model repeatedly. ``fjtag`` substitutes the struct and then the
CALLER'S OWN analyzer runs on it to completion: there is no callback seam and no
second pass, so the reusable unit is the transform/lift PAIR, not a driver.

Mirrors MATLAB ``matlab/src/solvers/TR/solver_tr_fjtag_analyzer.m``. It is a
mixin rather than a free function for the same reason ``ForkJoinDriverMixin``
is: python has no single result container, so the lift has to read and write
whatever the calling solver's own result object is.
"""

import numpy as np


class FJTagTransformMixin:
    """Tag-augment before the solve, fold the sibling classes back after it."""

    def _fjtag_require_network(self, solver_label):
        """The transform rewrites routing, so a bare NetworkStruct cannot serve.

        Kept separate from the expand so each caller can refuse in the order it
        always did, ahead of its own solver-specific refusals.
        """
        if not hasattr(self.model, 'copy') or not hasattr(self.model, 'get_linked_routing_matrix'):
            raise RuntimeError(
                "Native fork-join %s requires a Network model (not a bare NetworkStruct)."
                % solver_label)

    def _fjtag_expand(self, sn):
        """Substitute the tag-augmented struct and retain what the lift needs.

        Returns the augmented struct. The caller assigns it to ``self._sn`` and
        solves it exactly as it would an ordinary one.
        """
        from ..io.model_adapter import ModelAdapter
        korig = int(sn.nclasses)
        _fjmodel, fjsn, fjclassmap = ModelAdapter.fjtag(self.model)
        self._fj_foldback = (fjclassmap, korig, sn)
        return fjsn

    def _fjtag_lift(self, r):
        """Fold the auxiliary sibling classes back onto the original ones.

        Mutates ``r`` in place (Q, U, R, T, C, X), restores ``self._sn`` to the
        pre-augmentation struct and returns the arrival-rate matrix AN in
        ORIGINAL coordinates. The caller publishes AN only if its own result
        container carries an arrival-rate field: SolverSSA does, SolverCTMC has
        no such field and derives it later.
        """
        from ..api.fjnative import sn_fj_foldback
        from ..api.sn.getters import sn_get_arvr_from_tput, sn_pn_avg_rates
        fjclassmap, korig, orig_sn = self._fj_foldback
        r.Q, r.U, r.R, r.T, r.C, r.X = sn_fj_foldback(
            r.Q, r.U, r.R, r.T, r.C, r.X, fjclassmap, korig)
        # A Place counts tokens, not firings: rescale before deriving arrival
        # rates, and do it against the ORIGINAL struct whose class indices the
        # folded matrices are now in.
        r.T, _, r.R = sn_pn_avg_rates(orig_sn, r.Q, r.T, None, r.R)
        AN = np.atleast_2d(np.asarray(sn_get_arvr_from_tput(orig_sn, r.T, None), dtype=float))
        # A Join sees B sibling arrivals per released job, so its response time
        # is QN/AN and not QN/TN (the JMT convention).
        RN = np.array(r.R, dtype=float, copy=True)
        for ist in range(RN.shape[0]):
            for rr in range(min(korig, RN.shape[1])):
                if ist < AN.shape[0] and rr < AN.shape[1] and AN[ist, rr] > 0:
                    RN[ist, rr] = r.Q[ist, rr] / AN[ist, rr]
        r.R = RN
        self._sn = orig_sn
        return AN
