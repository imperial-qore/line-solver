"""Solver-agnostic driver of the fork-join fixed point.

The MMT (Dobre-Niu-Casale) and Heidelberger-Trivedi transformations rewrite a
fork-join network as a plain one: every fork becomes a router, every join a
zero-service delay, and the parallelism is carried by auxiliary open classes of
arrival rate ``(fanout-1)*forkLambda``. Each pass of the fixed point solves that
network, recomputes the synchronisation delays from the resulting metrics and
updates ``forkLambda``; the loop ends when the queue lengths stop moving.

Nothing in the loop reads a solver internal: the inner solve is reached only
through :meth:`ForkJoinDriverMixin._fj_inner_solver`, so any NetworkSolver that
can solve the transformed (fork-free) model can drive it. SolverMVA and SolverNC
are the two consumers, mirroring MATLAB ``@NetworkSolver/fjFixedPoint.m`` with
``@SolverMVA/mvaDispatch.m`` / ``@SolverNC/ncDispatch.m`` and the JAR
``jline.solvers.fj.FJFixedPoint`` with ``MVARunner.dispatch`` /
``SolverNC.ncDispatch``.
"""

import math
import os
import sys
import time
import warnings
from itertools import combinations
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from ..api.io.logging import line_debug, line_warning
from ..api.sn.network_struct import NodeType
from ..constants import GlobalConstants


class ForkJoinDriverMixin:
    """Fork-join transformation and fixed point, shared by MVA and NC."""

    def reset_fork_warm_start(self):
        """Discard the retained fork-join (MMT) fixed point."""
        self._fj_fork_lambda = None

    def _has_fork_join(self) -> bool:
        """Check if the model contains fork-join nodes."""
        if hasattr(self.model, 'has_fork') and callable(self.model.has_fork):
            return self.model.has_fork()
        # Fallback: check node types
        if self._sn is not None:
            from ..api.sn.network_struct import NodeType
            for nt in self._sn.nodetype:
                nt_val = nt.value if hasattr(nt, 'value') else int(nt)
                if nt_val == NodeType.FORK:
                    return True
        return False

    def _fj_network_type(self):
        """'open', 'closed' or 'mixed' for the model being driven.

        SolverMVA computes this into self.network_type when it extracts its
        parameters; SolverNC does not, so fall back to the struct populations.
        """
        nt = getattr(self, 'network_type', None)
        if nt is not None:
            return nt
        njobs = np.asarray(self._sn.njobs, dtype=float).ravel()
        has_open = bool(np.any(np.isinf(njobs)))
        has_closed = bool(np.any(np.isfinite(njobs)))
        if has_open and has_closed:
            return 'mixed'
        return 'closed' if has_closed else 'open'

    def _fj_publish(self, result):
        """Store the fork-join result in the solver's own result container.

        SolverMVA keeps a plain dict; SolverNC keeps a SolverNCReturn, and its
        getters read that, so the dict must be converted rather than assigned.
        """
        return result

    @staticmethod
    def _fj_result_dict(solver):
        """Inner-solve metrics as a plain dict, whatever container the solver uses.

        This is the python form of the callback contract that MATLAB
        (@SolverMVA/mvaDispatch.m) and the JAR (FJFixedPoint.InnerSolve) state
        explicitly: QN, UN, RN, TN, CN, XN, iter. SolverMVA already stores a
        dict; SolverNC stores a SolverNCReturn dataclass, whose fields carry the
        single-letter names.
        """
        res = getattr(solver, '_result', None)
        if res is None:
            return None
        if isinstance(res, dict):
            return res
        alias = {'QN': 'Q', 'UN': 'U', 'RN': 'R', 'TN': 'T', 'XN': 'X',
                 'CN': 'C', 'iter': 'it', 'lG': 'lG'}
        out = {}
        for key, attr in alias.items():
            val = getattr(res, attr, None)
            if val is None:
                val = getattr(res, key, None)
            if val is not None:
                out[key] = val
        return out

    def _fj_avg_tuple(self, solver):
        """(QN, UN, RN, TN, AN, WN) from an inner solve, container-independent."""
        avg = solver.getAvg()
        if avg is None:
            return None
        if isinstance(avg, (tuple, list)):
            return tuple(avg)
        res = self._fj_result_dict(solver)
        if res is None:
            return None
        shape = np.asarray(res['QN']).shape
        zeros = np.zeros(shape)
        return (res['QN'], res.get('UN', zeros), res.get('RN', zeros),
                res.get('TN', zeros), res.get('AN', zeros), res.get('WN', zeros))

    def _fj_inner_solver(self, nonfjmodel, method=None):
        """Build the solver applied to the transformed (fork-free) model.

        The default is SolverMVA, which is what the fixed point used when it
        lived inside SolverMVA. SolverNC overrides it so the same loop runs on
        the normalizing-constant analyzer.
        """
        from .solver_mva.solver_mva import SolverMVA
        if method is None:
            return SolverMVA(nonfjmodel)
        return SolverMVA(nonfjmodel, method=method)

    def _run_fork_join_analysis(self):
        """
        Run MVA analysis for fork-join networks.

        Uses the method specified in options.fork_join:
        - 'default' or 'mmt': MMT method (matches MATLAB default)
        - 'ht': Heidelberger-Trivedi method

        For closed fork-join networks:
        1. Transforms model to create auxiliary classes
        2. Solves transformed model with standard MVA
        3. Iteratively updates sync delays until convergence
        4. Merges auxiliary class results back into original classes

        For open fork-join networks:
        Returns None to fall back to standard analysis which then gets post-processed.
        """
        from ..api.sn.network_struct import NodeType
        from itertools import combinations

        sn = self._sn
        if sn is None:
            return None

        # MATLAB runAnalyzer.m line 68: checks self.model.hasFork, NOT sn.fj.
        # Fork-without-Join models (e.g., fj_nojoin) still need MMT transformation
        # to correctly compute throughputs with fork fanout.
        has_fork = any(int(nt) == int(NodeType.FORK) for nt in sn.nodetype)
        if not has_fork:
            return None

        # Get fork-join method from options
        # Default is 'mmt' to match MATLAB's default behavior
        fork_join_method = getattr(self.options, 'fork_join', 'default')
        if fork_join_method in ('default', 'mmt', 'fjt'):
            # Use MMT method (MATLAB default)
            return self._run_mmt_transformation()
        elif fork_join_method in ('ht', 'heidelberger-trivedi'):
            # Use Heidelberger-Trivedi method
            if self._fj_network_type() != 'closed':
                # H-T only works for closed networks
                return None
            return self._run_ht_transformation()
        else:
            # Unknown method, fall back to MMT
            return self._run_mmt_transformation()

    def _run_mmt_transformation(self):
        """
        Run MMT transformation for fork-join networks (MATLAB default method).

        For closed fork-join networks, MMT uses the H-T transformation but with
        a different sync delay formula: d0*fanOut - mean(ri) instead of d0*fanOut.

        This provides results consistent with MATLAB's default fork-join handling.

        Returns:
            Result dictionary or None if transformation fails
        """
        from ..io.model_adapter import ModelAdapter
        from ..api.sn.network_struct import NodeType
        from ..distributions.continuous import Exp
        from ..constants import GlobalConstants
        from itertools import combinations

        sn = self._sn
        if sn is None:
            return None

        # Check for Fork nodes (not sn.fj which only maps Fork→Join pairs)
        # Fork-without-Join models still need MMT transformation
        has_fork = any(int(nt) == int(NodeType.FORK) for nt in sn.nodetype)
        if not has_fork:
            return None

        # For closed networks, use H-T transformation with MMT sync delay formula
        if self._fj_network_type() == 'closed':
            return self._run_mmt_closed()
        else:
            # For open networks, MMT creates auxiliary open classes
            return self._run_mmt_open()

    def _run_mmt_closed(self):
        """
        Run MMT (Modified Method T) for closed fork-join networks.

        Matches MATLAB runAnalyzer.m lines 60-336 exactly:
        1. While forkLoop:
           a. If iter==1: call mmt() and sortForks()
           b. Else: update arrival rates, refreshRates()
           c. sn = nonfjmodel.getStruct(false)
           d. Check convergence on MERGED QN (from previous iteration)
           e. Run MVA on transformed model
           f. Fork handling: compute TNfork, update TN at join, compute sync delays
           g. Merge: save TN_orig, merge aux into original, restore TN, delete aux cols
        2. Return merged results
        """
        from ..io.model_adapter import ModelAdapter
        from ..api.sn.network_struct import NodeType
        from ..distributions.continuous import Exp
        from ..constants import GlobalConstants
        from ..api.sn.transforms import sn_get_residt_from_respt
        from ..api.sn.getters import sn_get_arvr_from_tput
        from itertools import combinations

        sn = self._sn
        orig_nclasses = sn.nclasses
        orig_nstations = sn.nstations

        # Initialize forkLambda for auxiliary classes (MATLAB line 63). Under an
        # outer iteration (e.g. SolverLN) the MMT fixed point is re-solved once
        # per outer iteration on a model whose parameters move only slightly;
        # restarting from FineTol each time discards the previous converged point
        # and makes the inner loop re-converge from scratch. Warm-start from the
        # retained iterate when it is conformant.
        num_forks = sum(1 for nt in sn.nodetype if int(nt) == int(NodeType.FORK))
        fork_lambda = GlobalConstants.FineTol * np.ones(2 * orig_nclasses * max(num_forks, 1))
        if getattr(self.options, 'fj_warmstart', True) \
                and getattr(self, '_fj_fork_lambda', None) is not None \
                and self._fj_fork_lambda.shape == fork_lambda.shape:
            fork_lambda = self._fj_fork_lambda.copy()

        # Initialize QN and QN_1 for convergence check (MATLAB lines 64-65)
        QN = GlobalConstants.Immediate * np.ones((1, orig_nclasses))
        QN_1 = np.zeros_like(QN)
        UN = np.zeros_like(QN)

        max_iter = getattr(self.options, 'max_iter', 1000)
        # The SOLVER's iteration tolerance, not GlobalConstants.CoarseTol.
        # CoarseTol is 1e-3 and is a display tolerance; using it here stops the
        # MMT iterate ~1e-3 short of its fixed point, so the auxiliary open
        # classes are merged back at an arrival rate that has not converged.
        # The symptom is a violation of flow conservation: on lqn_workflows the
        # activities on the branches of an AND fork come out 4.6e-4 above the
        # fork head and the join target, though each runs exactly once per
        # completion. The residual scales with this tolerance.
        tol = getattr(self.options, 'iter_tol', 1e-6)
        fork_loop = True
        fork_iter = 0
        total_iter = 0
        line_debug("Fork-join method: mmt", options=self.options)

        nonfjmodel = None
        fjclassmap = None
        fjforkmap = None
        fanout = None
        outer_forks = None
        parent_forks = None

        import os
        debug_mmt = os.environ.get('DEBUG_MMT', '0') == '1'

        while fork_loop and fork_iter < max_iter:
            fork_iter += 1
            line_debug("Fork-join iteration %d", fork_iter, options=self.options)

            # MATLAB lines 70-94: First iteration: transform model; else: update arrival rates
            if fork_iter == 1:
                # Reuse the cached MMT result when it was built from the very
                # same compiled struct, which is the case exactly when nothing
                # structural changed since (refresh_struct installs a new
                # NetworkStruct, so identity is the signal; holding the
                # reference also stops id() being recycled underneath us). The
                # fork topology the transformation encodes is a function of that
                # struct, so this is the precise condition for reuse. Rates may
                # well have changed, hence the refresh below.
                cache = getattr(self, '_mmt_cache', None)
                cache_usable = (
                    cache is not None
                    and cache.get('base_sn') is not None
                    and cache.get('base_sn') is getattr(self.model, '_sn', None)
                    and ModelAdapter.refresh_services_from_base(cache['mmt_result'])
                )
                if cache_usable:
                    mmt_result = cache['mmt_result']
                    nonfjmodel = mmt_result.nonfjmodel
                    fjclassmap = mmt_result.fjclassmap
                    fjforkmap = mmt_result.fjforkmap
                    fanout = mmt_result.fanout
                    outer_forks = cache['outer_forks']
                    parent_forks = cache['parent_forks']
                    # Only refresh rates on cached model (not full struct rebuild)
                    nonfjmodel.refresh_rates()
                else:
                    try:
                        mmt_result = ModelAdapter.mmt(self.model, fork_lambda[:orig_nclasses])
                    except Exception as e:
                        import warnings
                        warnings.warn(f"MMT transformation failed: {e}")
                        return None

                    nonfjmodel = mmt_result.nonfjmodel
                    fjclassmap = mmt_result.fjclassmap
                    fjforkmap = mmt_result.fjforkmap
                    fanout = mmt_result.fanout

                    if nonfjmodel is None or len(fjclassmap) == 0:
                        return None

                    # Sort forks using ModelAdapter (MATLAB line 76)
                    # Pass pre-computed routing matrix to avoid expensive refresh_struct
                    outer_forks, parent_forks = ModelAdapter.sort_forks(
                        sn, fjforkmap, fjclassmap, nonfjmodel,
                        routing_matrix=mmt_result.routing_matrix
                    )

                    # Cache for reuse across LN iterations, tagged with the
                    # struct it was derived from so a later structural rebuild
                    # cannot silently reuse a transformation of the old topology.
                    self._mmt_cache = {
                        'mmt_result': mmt_result,
                        'outer_forks': outer_forks,
                        'parent_forks': parent_forks,
                        'base_sn': getattr(self.model, '_sn', None),
                    }
            else:
                # Update arrival rates for auxiliary classes (MATLAB lines 78-94)
                source = nonfjmodel.get_source()
                if source is not None:
                    for r_aux in range(len(fjclassmap)):
                        s_orig = int(fjclassmap[r_aux])
                        if s_orig >= 0 and r_aux < len(fanout) and fanout[r_aux] > 0:
                            aux_class_idx = orig_nclasses + r_aux
                            if aux_class_idx < len(nonfjmodel._classes):
                                aux_class = nonfjmodel._classes[aux_class_idx]
                                arrival = source.get_arrival(aux_class)
                                if arrival is not None and not getattr(arrival, '_disabled', False):
                                    new_rate = (fanout[r_aux] - 1) * fork_lambda[r_aux]
                                    if new_rate > 0:
                                        source.set_arrival(aux_class, Exp(new_rate))
                # MATLAB uses refreshRates() here. NOTE: in practice this
                # usually does a FULL refresh_struct, because set_arrival above
                # nulls the cached struct. Preserving it across the mutation was
                # tried and measured as no help: the sync-delay block below still
                # rebuilds once per fork iteration (the AMVA solve invalidates
                # the struct), so refresh_struct stays ~97% of the fork loop
                # (dtmc_stochcomp -> scipy.linalg.inv). Fixing this needs the
                # struct preserved across the solve, not just across a mutation.
                nonfjmodel.refresh_rates()

            # MATLAB line 96: sn = nonfjmodel.getStruct(false)
            # Build struct if not yet available (iter 1 after relink clears it)
            if nonfjmodel._sn is None:
                nonfjmodel.refresh_struct()
            nonfjstruct = nonfjmodel._sn
            line_debug("Fork-join iter %d: rebuilt nonfjmodel struct (nstations=%d, nclasses=%d)",
                       fork_iter, nonfjstruct.nstations, nonfjstruct.nclasses, options=self.options)

            # MATLAB lines 97-108: convergence check on the full post-merge QN
            # from the previous iteration. MATLAB guards with
            # isequal(size(QN_1),size(QN)): a shape mismatch (first iterations,
            # before both tables have the augmented shape) is not converged.
            # Mixed absolute/relative test. A purely relative test cannot certify
            # convergence to zero: when the auxiliary class throughput at the join
            # is ~0, fork_lambda decays geometrically under the 0.5 damping, so QN
            # halves every iteration and the relative change stays pinned at 100%
            # down to denormals, running the loop to max_iter. An entry below the
            # absolute floor is numerically zero and has converged.
            if QN_1.shape != QN.shape:
                converged = False
            else:
                converged = bool(np.all(np.abs(QN_1 - QN) <=
                                        GlobalConstants.Zero + tol * np.abs(QN)))
            if debug_mmt and fork_iter <= 10:
                print(f"  [conv] iter={fork_iter}, converged={converged}, QN={QN.flatten()[:4]}, QN_1={QN_1.flatten()[:4]}", flush=True)
            if converged and fork_iter > 2:
                line_debug("Fork-join iter %d: converged (mixed abs/rel test)", fork_iter, options=self.options)
                fork_loop = False
                if debug_mmt:
                    print(f"Converged at iteration {fork_iter}")
            else:
                QN_1 = QN.copy()

            # Run MVA on the MMT-transformed model. MATLAB SolverMVA/runAnalyzer
            # forces the fork-join inner sub-solve to the AMVA linearizer
            # (method='amva') rather than exact BCMP mixed MVA (pfqn_mvamx); the
            # fork-join approximation is calibrated to the AMVA fixed point, and
            # exact mixed MVA gives a slightly different (3-6%) answer. Match MATLAB.
            nonfj_solver = self._fj_inner_solver(nonfjmodel, method='amva')
            nonfj_solver._skip_fork_join = True

            try:
                nonfj_result = self._fj_avg_tuple(nonfj_solver)
                if nonfj_result is None:
                    return None
                QN_new, UN_new, RN, TN, AN, WN = nonfj_result
                XN = nonfj_solver.getAvgSysTput()
            except Exception as e:
                import warnings
                warnings.warn(f"MVA failed: {e}")
                return None

            if debug_mmt and (fork_iter <= 5 or fork_iter % 50 == 0):
                print(f"\n=== MMT Iteration {fork_iter} ===")
                print(f"QN sum orig: {np.sum(QN_new[:, :orig_nclasses]):.4f}")

            # MATLAB lines 194-273: Fork handling - compute sync delays
            line_debug("Fork-join post-processing: method=%s, computing sync delays", "mmt", options=self.options)
            fork_indices = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.FORK)]

            # Use pre-computed routing matrix from MMT (avoids expensive refresh_struct)
            P_precomputed = mmt_result.routing_matrix

            for f in fork_indices:
                # MATLAB lines 200-206: Compute TNfork
                TNfork = np.zeros(orig_nclasses)
                for c in range(sn.nchains):
                    inchain = np.where(sn.chains[c, :] > 0)[0]
                    for r in inchain:
                        if parent_forks is not None and f < len(parent_forks):
                            pf = int(parent_forks[f])
                            if hasattr(sn, 'nodevisits') and sn.nodevisits is not None:
                                if c < len(sn.nodevisits) and sn.nodevisits[c] is not None:
                                    nv = sn.nodevisits[c]
                                    if pf < nv.shape[0] and r < nv.shape[1]:
                                        refstat = int(sn.refstat[r]) if r < len(sn.refstat) else 0
                                        visits_sum = 0
                                        if hasattr(sn, 'visits') and sn.visits is not None and c < len(sn.visits):
                                            v = sn.visits[c]
                                            stateful_idx = int(sn.stationToStateful[refstat]) if refstat < len(sn.stationToStateful) else refstat
                                            if v is not None and stateful_idx < v.shape[0]:
                                                visits_sum = np.sum(v[stateful_idx, inchain])
                                        if visits_sum > 0 and refstat < TN.shape[0]:
                                            TNfork[r] = (nv[pf, r] / visits_sum) * np.sum(TN[refstat, inchain])

                # MATLAB lines 208-209: Find join and auxiliary classes for this fork
                join_idx_arr = np.where(sn.fj[f, :] > 0)[0]
                if len(join_idx_arr) == 0:
                    continue
                join_idx = join_idx_arr[0]

                fork_aux_indices = np.where(fjforkmap == f)[0]

                for aux_idx in fork_aux_indices:
                    r = int(fjclassmap[aux_idx])
                    if r < 0:
                        continue

                    aux_class_idx = orig_nclasses + aux_idx

                    # MATLAB lines 212-217: Update forkLambda and TN at join
                    join_station = -1
                    if join_idx < len(nonfjstruct.nodeToStation):
                        join_station = int(nonfjstruct.nodeToStation[join_idx])

                    if join_station >= 0 and join_station < TN.shape[0] and r < TN.shape[1]:
                        # MATLAB line 215: TN(join,r) = TN(join,r) + sum(TN(join, find(fjclassmap==r))) - TN(join,s)
                        aux_indices_for_r = np.where(fjclassmap == r)[0]
                        aux_class_indices_for_r = orig_nclasses + aux_indices_for_r
                        tn_sum = TN[join_station, r] + np.sum(TN[join_station, aux_class_indices_for_r])
                        if aux_class_idx < TN.shape[1]:
                            tn_sum -= TN[join_station, aux_class_idx]
                        TN[join_station, r] = tn_sum
                        # MATLAB line 216: forkLambda(s) = mean([forkLambda(s); TN(join,r)])
                        fork_lambda[aux_idx] = np.mean([fork_lambda[aux_idx], tn_sum])

                    # MATLAB lines 218-220: Check if outer fork for class r
                    if f < outer_forks.shape[0] and r < outer_forks.shape[1]:
                        if outer_forks[f, r] == 0:
                            continue

                    # MATLAB lines 222-225: Find paths using findPathsCS
                    try:
                        P = P_precomputed
                        if P is not None:
                            P_combined = self._get_combined_routing_matrix(P, nonfjstruct, nonfjmodel)
                            to_merge = [r, aux_class_idx]
                            ri = ModelAdapter.find_paths_cs(
                                sn, P_combined, f, join_idx, r,
                                to_merge, QN_new, TN, 0.0,
                                fjclassmap, fjforkmap, nonfjmodel
                            )
                        else:
                            ri = np.array([GlobalConstants.FineTol, GlobalConstants.FineTol])
                    except Exception as e:
                        if debug_mmt:
                            print(f"    find_paths_cs exception: {e}")
                            import traceback
                            traceback.print_exc()
                        ri = np.array([GlobalConstants.FineTol, GlobalConstants.FineTol])

                    ri = np.maximum(ri, GlobalConstants.FineTol)

                    if len(ri) == 0:
                        continue

                    # MATLAB lines 226-232: Compute E[max] using inclusion-exclusion
                    lambdai = 1.0 / ri
                    d0 = 0.0
                    parallel_branches = len(ri)
                    for pow_val in range(parallel_branches):
                        combos = list(combinations(lambdai, pow_val + 1))
                        current_sum = sum(1.0 / sum(combo) for combo in combos)
                        d0 += ((-1) ** pow_val) * current_sum

                    # Get fanout (tasksPerLink) - MATLAB: sn.nodeparam{f}.fanOut
                    f_fanout = self._get_fanout(sn, f)

                    # MATLAB lines 234-235: Set sync delay
                    mean_ri = np.mean(ri)
                    sync_delay = max(d0 * f_fanout - mean_ri, GlobalConstants.FineTol)

                    if debug_mmt and (fork_iter <= 5 or fork_iter % 50 == 0):
                        print(f"  Fork {f}, r={r}: ri={ri}, d0={d0:.6f}, sync_delay={sync_delay:.6f}")

                    # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
                    saved_sn = nonfjmodel._sn
                    if aux_class_idx < len(nonfjmodel._classes):
                        nonfjmodel._nodes[join_idx].set_service(
                            nonfjmodel._classes[aux_class_idx], Exp.fit_mean(sync_delay)
                        )
                    if r < len(nonfjmodel._classes):
                        nonfjmodel._nodes[join_idx].set_service(
                            nonfjmodel._classes[r], Exp.fit_mean(sync_delay)
                        )
                    # MATLAB line 236: nonfjmodel.refreshRates()
                    # Restore the struct that set_service nulled. Go through the
                    # PUBLIC refresh_rates: when saved_sn is None it correctly
                    # falls back to refresh_struct and keeps _rates_dirty in
                    # step. Calling the private _refresh_rates here instead
                    # (a3648685c) crashed on `_sn.rates` whenever saved_sn was
                    # None. Measured: saved_sn is None on the first mutation of
                    # each fork iteration (the AMVA solve invalidates the
                    # struct), so this still rebuilds once per iteration.
                    nonfjmodel._sn = saved_sn
                    nonfjmodel.refresh_rates()

            # MATLAB lines 303-328: Merge auxiliary classes into original classes
            # This happens INSIDE the while loop, so convergence is checked on merged QN
            # Save TN at Join and Source stations (MATLAB line 304)
            join_indices_all = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.JOIN)]
            source_indices_all = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.SOURCE)]

            preserve_stations = []
            for j_idx in join_indices_all:
                if j_idx < len(nonfjstruct.nodeToStation):
                    preserve_stations.append(int(nonfjstruct.nodeToStation[j_idx]))
            for s_idx in source_indices_all:
                if s_idx < len(nonfjstruct.nodeToStation):
                    preserve_stations.append(int(nonfjstruct.nodeToStation[s_idx]))

            orig_class_set = sorted(set(int(c) for c in fjclassmap if c >= 0))

            TN_orig = {}
            for ist in preserve_stations:
                for oc in orig_class_set:
                    if ist < TN.shape[0] and oc < TN.shape[1]:
                        TN_orig[(ist, oc)] = TN[ist, oc]

            # MATLAB lines 306-319: merge back artificial classes
            for r_aux in range(len(fjclassmap)):
                s_orig = int(fjclassmap[r_aux])
                if s_orig >= 0:
                    aux_col = orig_nclasses + r_aux
                    if aux_col < QN_new.shape[1] and s_orig < QN_new.shape[1]:
                        QN_new[:, s_orig] += QN_new[:, aux_col]
                        UN_new[:, s_orig] += UN_new[:, aux_col]
                        TN[:, s_orig] += TN[:, aux_col]
                        with np.errstate(divide='ignore', invalid='ignore'):
                            RN[:, s_orig] = np.where(TN[:, s_orig] > 0,
                                                      QN_new[:, s_orig] / TN[:, s_orig], 0)

            # MATLAB line 321: Restore TN at Join and Source
            for (ist, oc), val in TN_orig.items():
                if ist < TN.shape[0] and oc < TN.shape[1]:
                    TN[ist, oc] = val

            # Update QN for the convergence check: post-merge originals only.
            # MATLAB deletes the auxiliary columns inside the loop
            # (QN(:,fjclassmap>0) = []), so its convergence table is exactly
            # the merged original-class block.
            QN = QN_new[:, :orig_nclasses].copy()
            UN = UN_new[:, :orig_nclasses].copy()

            _inner = self._fj_result_dict(nonfj_solver)
            total_iter += _inner.get('iter', 0) if _inner else 0
            if total_iter > 10000:
                break

        # The loop previously exhausted max_iter silently: convergence was only
        # ever reported through line_debug, so a non-converged MMT fixed point
        # was returned as a normal result.
        if fork_loop and fork_iter >= max_iter:
            line_warning('solver_mva',
                         'The fork-join (mmt) fixed point did not converge in options.iter_max=%d iterations; returning the interim solution.\n' % max_iter)
        # Retain the MMT iterate so that a subsequent runAnalyzer call on this
        # solver (an outer LN iteration) resumes the fixed point from here.
        self._fj_fork_lambda = np.asarray(fork_lambda).copy()

        # Extract final results: delete auxiliary class columns (MATLAB lines 323-328)
        QN_final = QN_new[:orig_nstations, :orig_nclasses].copy()
        UN_final = UN_new[:orig_nstations, :orig_nclasses].copy()
        RN_final = RN[:orig_nstations, :orig_nclasses].copy()
        TN_final = TN[:orig_nstations, :orig_nclasses].copy()
        XN_final = XN[:orig_nclasses].copy() if len(XN) >= orig_nclasses else XN.copy()

        # Zero out UN at Join stations
        for j_idx in join_indices_all:
            if j_idx < len(nonfjstruct.nodeToStation):
                join_ist = int(nonfjstruct.nodeToStation[j_idx])
                if join_ist < UN_final.shape[0]:
                    UN_final[join_ist, :] = 0

        # Compute arrival rates and residence times (MATLAB lines 341-342)
        AN_final = sn_get_arvr_from_tput(self._sn, TN_final, None)
        WN_final = sn_get_residt_from_respt(self._sn, RN_final, None)

        self._result = self._fj_publish({
            'QN': QN_final,
            'UN': UN_final,
            'RN': RN_final,
            'TN': TN_final,
            'AN': AN_final,
            'XN': XN_final,
            'WN': WN_final,
            'CN': np.sum(RN_final, axis=0),
            'runtime': 0.0,
            'method': 'mmt',
            'iter': total_iter
        })

        return self._result

    def _sort_forks(self, sn, fjforkmap, fjclassmap, nonfjmodel):
        """
        Sort forks to identify outer forks and their parents.
        Delegates to ModelAdapter.sort_forks which implements the full
        recursive nested fork detection matching MATLAB's sortForks.m.
        """
        from ..io.model_adapter import ModelAdapter
        return ModelAdapter.sort_forks(sn, fjforkmap, fjclassmap, nonfjmodel)

    @staticmethod
    def _get_fanout(sn, f):
        """Get fanOut (tasksPerLink) for fork node f from sn.nodeparam."""
        f_fanout = 1  # Default: 1 task per link
        if hasattr(sn, 'nodeparam') and sn.nodeparam is not None:
            fp = None
            if isinstance(sn.nodeparam, dict) and f in sn.nodeparam:
                fp = sn.nodeparam[f]
            elif isinstance(sn.nodeparam, list) and f < len(sn.nodeparam):
                fp = sn.nodeparam[f]
            if fp is not None:
                if isinstance(fp, dict) and 'fanOut' in fp:
                    f_fanout = int(fp['fanOut'])
                elif hasattr(fp, 'fanOut'):
                    f_fanout = int(fp.fanOut)
        return f_fanout

    def _get_combined_routing_matrix(self, P, sn, nonfjmodel=None):
        """
        Convert cell-based routing matrix P to combined format for find_paths_cs.

        The combined matrix has dimensions (nclasses * nnodes) x (nclasses * nnodes).
        Entry at (r*nnodes + i, s*nnodes + j) represents routing probability
        from (class r, node i) to (class s, node j).

        ClassSwitch nodes are resolved: their class-switching behavior is folded
        into the routing probabilities so that routing bypasses CS nodes. This
        matches MATLAB's behavior where cell2mat(getLinkedRoutingMatrix) produces
        a combined matrix with CS nodes already resolved.

        Args:
            P: Cell-based routing matrix P[r][s] = nnodes x nnodes matrix
            sn: Network structure (for extracting node count)
            nonfjmodel: Optional transformed model for ClassSwitch resolution

        Returns:
            Combined routing matrix as numpy array
        """
        from ..lang.nodes import ClassSwitch as CSNode

        if P is None:
            return None

        nclasses = len(P)
        nnodes = P[0][0].shape[0] if nclasses > 0 and len(P[0]) > 0 else 0

        # Create combined matrix of size (nclasses*nnodes, nclasses*nnodes)
        P_combined = np.zeros((nclasses * nnodes, nclasses * nnodes))
        for r in range(nclasses):
            for s in range(nclasses):
                if r < len(P) and s < len(P[r]):
                    P_combined[r*nnodes:(r+1)*nnodes, s*nnodes:(s+1)*nnodes] = P[r][s]

        # Resolve ClassSwitch nodes: fold class-switching into routing
        # probabilities so that find_paths_cs can traverse paths correctly.
        # In Python, the linked routing matrix routes (class_r, predecessor) ->
        # (class_r, CS_node) and then (switched_class, CS_node) -> (switched_class,
        # successor). MATLAB resolves this into direct (class_r, predecessor) ->
        # (switched_class, successor) routing.
        if nonfjmodel is None:
            nonfjmodel = self.model
        cs_nodes = []
        if nonfjmodel is not None and hasattr(nonfjmodel, '_nodes'):
            for node_idx, node in enumerate(nonfjmodel._nodes):
                if isinstance(node, CSNode) and node_idx < nnodes:
                    switch_matrix = getattr(node, '_switch_matrix', None)
                    if switch_matrix is not None:
                        cs_nodes.append((node_idx, switch_matrix))

        for cs_idx, switch_matrix in cs_nodes:
            # For each CS node, resolve: incoming * switch * outgoing
            for r_in in range(nclasses):
                # Incoming: P_combined[r_in*n + i, r_in*n + cs] for all predecessor i
                in_col = r_in * nnodes + cs_idx
                incoming = P_combined[:, in_col].copy()
                if np.sum(incoming) == 0:
                    continue

                # Determine switched class using switch matrix
                for r_out in range(nclasses):
                    if r_out >= switch_matrix.shape[0] or r_in >= switch_matrix.shape[1]:
                        continue
                    # switch_matrix[r_in, r_out] = prob of switching from r_in to r_out
                    sw_prob = switch_matrix[r_in, r_out]
                    if sw_prob <= 0:
                        continue

                    # Outgoing: P_combined[r_out*n + cs, r_out*n + j] for all successor j
                    out_row = r_out * nnodes + cs_idx
                    outgoing = P_combined[out_row, :].copy()
                    if np.sum(outgoing) == 0:
                        continue

                    # Add resolved routing: predecessor -> successor with class switch
                    for pred_idx in np.where(incoming > 0)[0]:
                        for succ_idx in np.where(outgoing > 0)[0]:
                            P_combined[pred_idx, succ_idx] += (
                                incoming[pred_idx] * sw_prob * outgoing[succ_idx]
                            )

            # Zero out CS node rows and columns in the combined matrix
            for r in range(nclasses):
                P_combined[r * nnodes + cs_idx, :] = 0
                P_combined[:, r * nnodes + cs_idx] = 0

        return P_combined

    def _merge_mmt_results(self, QN, UN, RN, TN, XN, orig_sn, nonfjstruct, fjclassmap, fjforkmap, orig_nclasses, orig_nstations):
        """
        Merge results from MMT transformed model back to original model structure.
        """
        from ..api.sn.network_struct import NodeType
        from ..api.sn.transforms import sn_get_residt_from_respt

        # Find join station indices
        join_stations = []
        fork_indices = [i for i in range(orig_sn.nnodes) if int(orig_sn.nodetype[i]) == int(NodeType.FORK)]
        for f in fork_indices:
            join_idx_arr = np.where(orig_sn.fj[f, :] > 0)[0]
            if len(join_idx_arr) > 0:
                join_idx = join_idx_arr[0]
                if join_idx < len(nonfjstruct.nodeToStation):
                    join_stations.append(int(nonfjstruct.nodeToStation[join_idx]))

        # Find Source station index
        source_stations = []
        source_indices = [i for i in range(orig_sn.nnodes) if int(orig_sn.nodetype[i]) == int(NodeType.SOURCE)]
        for src_idx in source_indices:
            if src_idx < len(nonfjstruct.nodeToStation):
                source_stations.append(int(nonfjstruct.nodeToStation[src_idx]))

        # Merge auxiliary class metrics into original classes
        QN_merged = QN.copy()
        UN_merged = UN.copy()
        RN_merged = RN.copy()
        TN_merged = TN.copy()

        # Save TN at Join and Source stations before merging (MATLAB lines 304, 321)
        # These throughputs should NOT be summed - restore them after merge
        special_stations = join_stations + source_stations
        orig_class_indices = np.unique([int(fjclassmap[i]) for i in range(len(fjclassmap)) if fjclassmap[i] >= 0])
        TN_orig = {}
        for st in special_stations:
            if st < TN.shape[0]:
                for oc in orig_class_indices:
                    if oc < TN.shape[1]:
                        TN_orig[(st, oc)] = TN[st, oc]

        # Merge auxiliary class metrics into original classes
        # Auxiliary classes are at indices orig_nclasses + aux_idx in the transformed model
        # MATLAB runAnalyzer.m lines 306-316: merge and recompute RN = QN / TN
        for aux_idx in range(len(fjclassmap)):
            orig_idx = fjclassmap[aux_idx]
            aux_class_idx = orig_nclasses + aux_idx  # Actual index of auxiliary class
            if orig_idx >= 0 and aux_class_idx < QN_merged.shape[1] and orig_idx < QN_merged.shape[1]:
                QN_merged[:, orig_idx] += QN_merged[:, aux_class_idx]
                UN_merged[:, orig_idx] += UN_merged[:, aux_class_idx]
                TN_merged[:, orig_idx] += TN_merged[:, aux_class_idx]
                # Recompute RN = QN / TN after merging (MATLAB line 316)
                with np.errstate(divide='ignore', invalid='ignore'):
                    RN_merged[:, orig_idx] = np.where(TN_merged[:, orig_idx] > 0,
                                                       QN_merged[:, orig_idx] / TN_merged[:, orig_idx], 0)

        # Restore TN at Join and Source stations (MATLAB line 321)
        # These throughputs should NOT be summed - restore them after merge
        for st, oc in TN_orig.keys():
            if st < TN_merged.shape[0] and oc < TN_merged.shape[1]:
                TN_merged[st, oc] = TN_orig[(st, oc)]

        # Zero out UN at Join stations
        for js in join_stations:
            if js < UN_merged.shape[0]:
                UN_merged[js, :orig_nclasses] = 0

        # Extract only original classes and original stations
        # Remove auxiliary delay stations if any
        stations_to_keep = list(range(min(orig_nstations, QN_merged.shape[0])))

        QN_final = QN_merged[stations_to_keep, :][:, :orig_nclasses]
        UN_final = UN_merged[stations_to_keep, :][:, :orig_nclasses]
        RN_final = RN_merged[stations_to_keep, :][:, :orig_nclasses]
        TN_final = TN_merged[stations_to_keep, :][:, :orig_nclasses]
        XN_final = XN[:orig_nclasses].copy() if len(XN) >= orig_nclasses else XN.copy()

        # Compute arrival rates using sn_get_arvr_from_tput (MATLAB runAnalyzer.m line 341)
        from ..api.sn.getters import sn_get_arvr_from_tput
        AN_final = sn_get_arvr_from_tput(self._sn, TN_final, None)

        # Compute residence times
        WN_final = sn_get_residt_from_respt(self._sn, RN_final, None)

        # Store results
        self._result = self._fj_publish({
            'QN': QN_final,
            'UN': UN_final,
            'RN': RN_final,
            'TN': TN_final,
            'AN': AN_final,
            'XN': XN_final,
            'WN': WN_final,
            'CN': np.sum(RN_final, axis=0),
            'runtime': 0.0,
            'method': 'mmt',
            'iter': iter_count if 'iter_count' in dir() else 0
        })

        return self._result

    def _merge_ht_results(self, QN, UN, RN, TN, XN, orig_sn, nonfjstruct, fjclassmap, fj_auxiliary_delays, orig_nclasses, orig_nstations):
        """
        Merge results from H-T/MMT transformed model back to original model structure.

        This function:
        1. Finds join and auxiliary delay station indices
        2. Zeros out original class metrics at join stations
        3. Removes auxiliary delay station rows
        4. Merges auxiliary class metrics into original classes
        5. Computes final metrics and stores result
        """
        from ..api.sn.network_struct import NodeType
        from ..api.sn.transforms import sn_get_residt_from_respt

        # Find fork indices in original model
        fork_indices = []
        for i, nt in enumerate(orig_sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK.value if hasattr(NodeType.FORK, 'value') else NodeType.FORK:
                fork_indices.append(i)

        # Find join and aux delay station indices
        join_stations = []
        aux_delay_stations = []
        for f in fork_indices:
            join_idx_arr = np.where(orig_sn.fj[f, :] > 0)[0]
            if len(join_idx_arr) > 0:
                join_idx = join_idx_arr[0]
                if join_idx < len(nonfjstruct.nodeToStation):
                    join_stations.append(nonfjstruct.nodeToStation[join_idx])
                aux_delay_idx = fj_auxiliary_delays.get(join_idx, None)
                if aux_delay_idx is not None and aux_delay_idx < len(nonfjstruct.nodeToStation):
                    aux_delay_stations.append(nonfjstruct.nodeToStation[aux_delay_idx])

        # Get original class indices that have auxiliary classes
        orig_classes_with_aux = sorted(set(fjclassmap[fjclassmap >= 0]))

        # Save original class throughputs at join stations before zeroing
        TN_orig_join = {}
        for js in join_stations:
            for oc in orig_classes_with_aux:
                if js < TN.shape[0] and oc < TN.shape[1]:
                    TN_orig_join[(js, oc)] = TN[js, oc]

        # Zero out original class metrics at join stations
        for js in join_stations:
            for oc in orig_classes_with_aux:
                if js < QN.shape[0] and oc < QN.shape[1]:
                    QN[js, oc] = 0
                    RN[js, oc] = 0
                    TN[js, oc] = 0
                    UN[js, oc] = 0

        # Remove auxiliary delay station rows (create mask for stations to keep)
        stations_to_keep = [i for i in range(QN.shape[0]) if i not in aux_delay_stations]

        QN_merged = QN[stations_to_keep, :]
        UN_merged = UN[stations_to_keep, :]
        RN_merged = RN[stations_to_keep, :]
        TN_merged = TN[stations_to_keep, :]

        # Merge auxiliary classes into original classes
        # fjclassmap[class_idx] = -1 for original classes, = orig_class_idx for auxiliary classes
        # Auxiliary classes are stored at their actual indices in the fjclassmap
        for aux_class_idx in range(len(fjclassmap)):
            orig_idx = fjclassmap[aux_class_idx]
            # Only process auxiliary classes (those that map to an original class)
            if orig_idx >= 0 and aux_class_idx < QN_merged.shape[1] and orig_idx < QN_merged.shape[1]:
                QN_merged[:, orig_idx] += QN_merged[:, aux_class_idx]
                UN_merged[:, orig_idx] += UN_merged[:, aux_class_idx]
                # Add all throughputs of the auxiliary classes to facilitate the computation of the response times
                TN_merged[:, orig_idx] += TN_merged[:, aux_class_idx]
                # Recompute RN = QN / TN after merging (matching JAR behavior)
                for i in range(RN_merged.shape[0]):
                    if TN_merged[i, orig_idx] > 0:
                        RN_merged[i, orig_idx] = QN_merged[i, orig_idx] / TN_merged[i, orig_idx]

        # Restore original class throughputs at join stations
        for (old_js, oc), tput in TN_orig_join.items():
            if old_js in stations_to_keep:
                new_js_idx = stations_to_keep.index(old_js)
                if new_js_idx < TN_merged.shape[0] and oc < TN_merged.shape[1]:
                    TN_merged[new_js_idx, oc] = tput

        # Zero out UN at Join stations for all original classes
        # Join has Immediate service in the original model, so utilization must be 0
        for js in join_stations:
            if js in stations_to_keep:
                new_js_idx = stations_to_keep.index(js)
                if new_js_idx < UN_merged.shape[0]:
                    UN_merged[new_js_idx, :orig_nclasses] = 0

        # Extract only original classes
        QN_final = QN_merged[:, :orig_nclasses]
        UN_final = UN_merged[:, :orig_nclasses]
        RN_final = RN_merged[:, :orig_nclasses]
        TN_final = TN_merged[:, :orig_nclasses]
        XN_final = XN[:orig_nclasses].copy()

        # Compute arrival rates using sn_get_arvr_from_tput (MATLAB runAnalyzer.m line 341)
        from ..api.sn.getters import sn_get_arvr_from_tput
        AN_final = sn_get_arvr_from_tput(self._sn, TN_final, None)

        # Compute residence times from response times (WN = RN * V)
        WN_final = sn_get_residt_from_respt(self._sn, RN_final, None)

        # Store results
        self._result = self._fj_publish({
            'QN': QN_final,
            'UN': UN_final,
            'RN': RN_final,
            'TN': TN_final,
            'AN': AN_final,
            'XN': XN_final,
            'WN': WN_final,
            'CN': np.sum(RN_final, axis=0),
            'runtime': 0.0,
            'method': 'mmt',
            'iter': 0
        })

        return self._result

    def _run_mmt_open(self):
        """
        Run MMT for open fork-join networks.

        1-to-1 port of MATLAB runAnalyzer.m MMT implementation:
        1. Transform model using ModelAdapter.mmt() to create auxiliary classes
        2. Iterate until convergence:
           - Update auxiliary class arrival rates based on fanout
           - Solve transformed model using standard MVA
           - Compute TNfork and update forkLambda
           - Call findPathsCS to get response times on parallel paths
           - Compute sync delay d0 and update Join service times
        3. Merge auxiliary class results back into original classes

        Returns:
            Result dictionary or None if transformation fails
        """
        from ..io.model_adapter import ModelAdapter
        from ..api.sn.network_struct import NodeType
        from ..api.sn.transforms import sn_get_residt_from_respt
        from ..distributions.continuous import Exp
        from ..constants import GlobalConstants
        from itertools import combinations

        sn = self._sn
        if sn is None:
            return None

        orig_nclasses = sn.nclasses
        orig_nstations = sn.nstations
        orig_sn = sn  # Save original sn

        # Find fork indices
        fork_indices = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.FORK)]
        num_forks = len(fork_indices)
        if num_forks == 0:
            return None

        # Initialize forkLambda (MATLAB line 63)
        fork_lambda = GlobalConstants.FineTol * np.ones(2 * orig_nclasses * num_forks)

        # Initialize QN and QN_1 for convergence check (MATLAB lines 64-65)
        QN = GlobalConstants.Immediate * np.ones(orig_nclasses)
        QN_1 = np.zeros(orig_nclasses)
        UN = np.zeros(orig_nclasses)

        fork_loop = True
        fork_iter = 0
        stop_next = False  # convergence latch: one more solve after firing
        max_iter = 500  # options.iter_max

        nonfjmodel = None
        fjclassmap = None
        fjforkmap = None
        fanout = None
        outer_forks = None
        parent_forks = None
        total_iter = 0

        DEBUG_MMT = os.environ.get("DEBUG_MMT", "0") == "1"  # TEMP DEBUG
        while fork_loop and fork_iter < max_iter:
            fork_iter += 1
            total_iter = fork_iter

            if DEBUG_MMT and fork_iter <= 3:
                print(f"DEBUG: Iteration {fork_iter}, fork_lambda[:5] = {fork_lambda[:5]}")

            # First iteration: transform model (MATLAB lines 70-77)
            if fork_iter == 1:
                try:
                    mmt_result = ModelAdapter.mmt(self.model, fork_lambda[:orig_nclasses])
                    nonfjmodel = mmt_result.nonfjmodel
                    fjclassmap = mmt_result.fjclassmap
                    fjforkmap = mmt_result.fjforkmap
                    fanout = mmt_result.fanout

                    if len(fjclassmap) == 0:
                        return None

                    # Sort forks (MATLAB line 76)
                    outer_forks, parent_forks = ModelAdapter.sort_forks(
                        sn, fjforkmap, fjclassmap, nonfjmodel
                    )
                except Exception as e:
                    import warnings
                    warnings.warn(f"MMT transformation failed: {e}. Falling back to standard MVA.")
                    return None
            else:
                # Subsequent iterations: update auxiliary class arrival rates (MATLAB lines 78-94)
                try:
                    from ..distributions.continuous import Disabled
                    nonfj_source = nonfjmodel.get_source()
                    nonfj_classes = nonfjmodel._classes
                    if DEBUG_MMT and fork_iter <= 3:
                        print(f"DEBUG: len(fjclassmap)={len(fjclassmap)}, fjclassmap={fjclassmap}, fanout={fanout}")
                        print(f"DEBUG: nonfj_classes count={len(nonfj_classes)}")
                    for r in range(len(fjclassmap)):  # r is auxiliary class index
                        s = int(fjclassmap[r])  # s is original class
                        if DEBUG_MMT and fork_iter <= 3:
                            print(f"DEBUG: r={r}, s={s}, check s>=0:{s>=0}, r<len(fanout):{r<len(fanout)}, fanout[r]>0:{fanout[r]>0 if r<len(fanout) else 'N/A'}")
                        if s >= 0 and r < len(fanout) and fanout[r] > 0:
                            # MATLAB: nonfjSource.arrivalProcess{r}.setRate((fanout(r)-1)*forkLambda(r))
                            # In Python, use set_arrival with the aux class object
                            aux_class_idx = orig_nclasses + r
                            if aux_class_idx < len(nonfj_classes):
                                aux_class = nonfj_classes[aux_class_idx]
                                arr_proc = nonfj_source.get_arrival(aux_class)
                                is_disabled = isinstance(arr_proc, Disabled) if arr_proc is not None else True
                                if DEBUG_MMT and fork_iter <= 3:
                                    print(f"DEBUG: aux_class={aux_class.name}, arr_proc={arr_proc}, is_disabled={is_disabled}")
                                if arr_proc is not None and not is_disabled:
                                    new_rate = (fanout[r] - 1) * fork_lambda[r]
                                    if DEBUG_MMT and fork_iter <= 3:
                                        print(f"DEBUG: Updating aux class {aux_class.name} rate to {new_rate}, fanout={fanout[r]}, forkLambda={fork_lambda[r]}")
                                    if new_rate > 0:
                                        # Update by setting a new Exp distribution with the updated rate
                                        nonfj_source.set_arrival(aux_class, Exp(new_rate))
                    # Refresh the model struct after updating arrival rates
                    nonfjmodel.refresh_struct()
                except Exception as ex:
                    if DEBUG_MMT:
                        print(f"DEBUG: Exception in arrival rate update: {ex}")
                        import traceback
                        traceback.print_exc()

            # Solve nonfjmodel (MATLAB line 96, then lines 181-190)
            nonfjmodel.refresh_struct()
            nonfj_sn = nonfjmodel._sn

            from . import SolverMVA
            transformed_solver = self._fj_inner_solver(nonfjmodel)
            transformed_solver._skip_fork_join = True  # Avoid recursion
            transformed_solver._force_method = 'amva'

            try:
                transformed_solver.runAnalyzer()
                trans_result = self._fj_result_dict(transformed_solver)
            except Exception as e:
                import warnings
                warnings.warn(f"MVA on transformed model failed: {e}.")
                break

            if trans_result is None:
                break

            # Get results from transformed model
            QN_full = trans_result['QN'].copy()
            UN_full = trans_result['UN'].copy()
            RN_full = trans_result['RN'].copy()
            TN_full = trans_result['TN'].copy()
            CN_full = trans_result.get('CN', np.zeros(QN_full.shape[1])).copy()
            XN_full = trans_result['XN'].copy()

            if DEBUG_MMT and fork_iter <= 3:
                print(f"DEBUG: After solve iter {fork_iter}, TN_full shape={TN_full.shape}")
                print(f"DEBUG: TN_full =\n{TN_full}")

            # MATLAB lines 97-108: convergence check, mirrored exactly.
            # MATLAB (a) folds the auxiliary columns into their original classes
            # in-loop and tests the (stations x Korig) table ELEMENTWISE
            # (0/0 -> converged entry, x/0 -> not converged), and (b) performs
            # the check at the TOP of the iteration, before the solve, so one
            # more solve runs after convergence fires: the endpoint is that
            # final solve's tables. Compare the merged tables of the two
            # PREVIOUS solves here, then always fall through to the solve.
            if stop_next:
                # Convergence fired on the previous pair of solves; this
                # iteration's solve is the final one (MATLAB's endpoint).
                fork_loop = False
            else:
                QN_check = QN_full[:, :orig_nclasses].copy()
                for a_idx in range(len(fjclassmap)):
                    s0 = int(fjclassmap[a_idx])
                    aux_col = orig_nclasses + a_idx
                    if s0 >= 0 and aux_col < QN_full.shape[1]:
                        QN_check[:, s0] += QN_full[:, aux_col]
                # Mixed absolute/relative test; see the closed-model loop above.
                # A relative-only test never certifies convergence to zero, so a
                # geometrically decaying auxiliary iterate ran to max_iter.
                if not (isinstance(QN_1, np.ndarray) and QN_1.shape == QN_check.shape):
                    qn_converged = False
                else:
                    qn_converged = bool(np.all(
                        np.abs(QN_1 - QN_check) <=
                        GlobalConstants.Zero +
                        getattr(self.options, 'iter_tol', 1e-6) * np.abs(QN_check)))
                if DEBUG_MMT:
                    print(f"OPENCONV {fork_iter} {qn_converged}", flush=True)
                if qn_converged and fork_iter > 2:
                    stop_next = True
                else:
                    QN_1 = QN_check.copy()

            # MATLAB lines 194-237: Compute sync delays for each fork
            # nonfjstruct = nonfj_sn (the transformed model's sn)
            nonfjstruct = nonfj_sn

            # Use pre-computed routing matrix from MMT (avoids expensive refresh_struct)
            P_precomputed = mmt_result.routing_matrix

            for f in fork_indices:
                # MATLAB lines 200-206: Compute TNfork
                TNfork = np.zeros(orig_nclasses)
                for c in range(sn.nchains):
                    inchain = np.where(sn.chains[c, :] > 0)[0]
                    for r in inchain:
                        if parent_forks is not None and f < len(parent_forks):
                            pf = int(parent_forks[f])
                            if hasattr(sn, 'nodevisits') and sn.nodevisits is not None:
                                if c < len(sn.nodevisits) and sn.nodevisits[c] is not None:
                                    nv = sn.nodevisits[c]
                                    if pf < nv.shape[0] and r < nv.shape[1]:
                                        visits_sum = 0
                                        refstat = int(sn.refstat[r]) if r < len(sn.refstat) else 0
                                        if hasattr(sn, 'visits') and sn.visits is not None and c < len(sn.visits):
                                            v = sn.visits[c]
                                            stateful_idx = int(sn.stationToStateful[refstat]) if refstat < len(sn.stationToStateful) else refstat
                                            if v is not None and stateful_idx < v.shape[0]:
                                                visits_sum = np.sum(v[stateful_idx, inchain])
                                        if visits_sum > 0 and refstat < TN_full.shape[0]:
                                            TNfork[r] = (nv[pf, r] / visits_sum) * np.sum(TN_full[refstat, inchain])

                # MATLAB lines 208-209: Find join and auxiliary classes for this fork
                join_idx_arr = np.where(sn.fj[f, :] > 0)[0]
                has_join = len(join_idx_arr) > 0
                join_idx = join_idx_arr[0] if has_join else None

                fork_aux_classes = np.where(fjforkmap == f)[0]

                for s_idx in fork_aux_classes:
                    r = int(fjclassmap[s_idx])  # original class
                    if r < 0 or r >= orig_nclasses:
                        continue

                    # MATLAB lines 212-217: Update forkLambda
                    if not has_join:
                        # No join: use TNfork (MATLAB line 213)
                        fork_lambda[s_idx] = 0.5 * (fork_lambda[s_idx] + TNfork[r])
                    else:
                        join_ist = int(sn.nodeToStation[join_idx]) if join_idx < len(sn.nodeToStation) else -1
                        if join_ist >= 0 and join_ist < TN_full.shape[0]:
                            # TN at join for class r plus auxiliary classes mapped to r
                            aux_for_r = np.where(fjclassmap == r)[0]
                            tn_join_r = TN_full[join_ist, r] + np.sum(TN_full[join_ist, orig_nclasses + aux_for_r]) - TN_full[join_ist, orig_nclasses + s_idx]
                            fork_lambda[s_idx] = 0.5 * (fork_lambda[s_idx] + tn_join_r)

                    # MATLAB line 218: if isempty(joinIdx) || ~outer_forks(f, r) -> skip sync delay
                    if DEBUG_MMT and fork_iter <= 3:
                        print(f"DEBUG sync: f={f}, r={r}, has_join={has_join}, outer_forks shape={outer_forks.shape if outer_forks is not None else None}")
                        if outer_forks is not None and f < outer_forks.shape[0] and r < outer_forks.shape[1]:
                            print(f"DEBUG sync: outer_forks[{f},{r}]={outer_forks[f, r]}")
                    if not has_join or outer_forks is None or (f < outer_forks.shape[0] and r < outer_forks.shape[1] and not outer_forks[f, r]):
                        if DEBUG_MMT and fork_iter <= 3:
                            print(f"DEBUG sync: SKIPPING sync delay for f={f}, r={r}")
                        continue

                    # MATLAB lines 222-225: Find paths using findPathsCS
                    try:
                        P = P_precomputed
                        if DEBUG_MMT and fork_iter <= 3:
                            print(f"DEBUG findPaths: P is None: {P is None}")
                        if P is not None:
                            # Convert cell array to combined matrix
                            P_combined = SolverMVA._combine_routing_matrix(P)
                            if DEBUG_MMT and fork_iter <= 3:
                                print(f"DEBUG findPaths: P_combined shape={P_combined.shape if P_combined is not None else None}")
                            ri = ModelAdapter.find_paths_cs(
                                sn, P_combined, f, join_idx, r,
                                [r, orig_nclasses + s_idx], QN_full, TN_full, 0,
                                fjclassmap, fjforkmap, nonfjmodel
                            )
                            if DEBUG_MMT and fork_iter <= 3:
                                print(f"DEBUG findPaths: ri={ri}, len(ri)={len(ri)}")
                        else:
                            ri = np.array([])
                    except Exception as ex:
                        if DEBUG_MMT and fork_iter <= 3:
                            print(f"DEBUG findPaths: Exception {ex}")
                            import traceback
                            traceback.print_exc()
                        ri = np.array([])

                    if len(ri) < 2:
                        if DEBUG_MMT and fork_iter <= 3:
                            print(f"DEBUG findPaths: SKIPPING sync delay, len(ri)={len(ri)} < 2")
                        continue

                    # MATLAB lines 226-232: Compute d0 using order statistics
                    lambdai = 1.0 / ri
                    d0 = 0.0
                    parallel_branches = len(ri)
                    for pow_k in range(parallel_branches):
                        current_sum = 0.0
                        for subset in combinations(range(parallel_branches), pow_k + 1):
                            lambda_sum = np.sum(lambdai[list(subset)])
                            if lambda_sum > 0:
                                current_sum += 1.0 / lambda_sum
                        d0 += ((-1) ** pow_k) * current_sum

                    # MATLAB lines 234-236: Set sync delay at Join
                    fan_out = SolverMVA._get_fanout(sn, f)

                    sync_delay = max(0, d0 * fan_out - np.mean(ri))
                    if DEBUG_MMT and fork_iter <= 3:
                        print(f"DEBUG sync_delay: d0={d0}, fan_out={fan_out}, mean(ri)={np.mean(ri)}, sync_delay={sync_delay}")

                    try:
                        nonfj_nodes = nonfjmodel._nodes
                        nonfj_classes = nonfjmodel._classes
                        if join_idx < len(nonfj_nodes):
                            join_node = nonfj_nodes[join_idx]
                            if DEBUG_MMT and fork_iter <= 3:
                                print(f"DEBUG sync_delay: join_idx={join_idx}, join_node type={type(join_node)}, has set_service={hasattr(join_node, 'set_service')}")
                            if sync_delay > 0 and hasattr(join_node, 'set_service'):
                                # Set for auxiliary class
                                if orig_nclasses + s_idx < len(nonfj_classes):
                                    if DEBUG_MMT and fork_iter <= 3:
                                        print(f"DEBUG sync_delay: Setting aux class {nonfj_classes[orig_nclasses + s_idx].name} sync delay={sync_delay}")
                                    join_node.set_service(nonfj_classes[orig_nclasses + s_idx], Exp.fit_mean(sync_delay))
                                # Set for original class
                                if r < len(nonfj_classes):
                                    if DEBUG_MMT and fork_iter <= 3:
                                        print(f"DEBUG sync_delay: Setting orig class {nonfj_classes[r].name} sync delay={sync_delay}")
                                    join_node.set_service(nonfj_classes[r], Exp.fit_mean(sync_delay))
                                nonfjmodel.refresh_struct()
                    except Exception as ex:
                        if DEBUG_MMT:
                            print(f"DEBUG sync_delay: Exception {ex}")
                            import traceback
                            traceback.print_exc()

        # After loop: merge results (MATLAB lines 303-328)
        QN = QN_full.copy()
        UN = UN_full.copy()
        RN = RN_full.copy()
        TN = TN_full.copy()
        CN = CN_full.copy()
        XN = XN_full.copy()

        # MATLAB line 304: Save TN_orig for Join and Source
        join_indices = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.JOIN)]
        source_indices = [i for i in range(sn.nnodes) if int(sn.nodetype[i]) == int(NodeType.SOURCE)]

        preserve_ist_list = []
        for j_idx in join_indices:
            if j_idx < len(nonfj_sn.nodeToStation):
                preserve_ist_list.append(int(nonfj_sn.nodeToStation[j_idx]))
        for s_idx in source_indices:
            if s_idx < len(nonfj_sn.nodeToStation):
                preserve_ist_list.append(int(nonfj_sn.nodeToStation[s_idx]))

        orig_classes_with_aux = [int(s) for s in set(fjclassmap) if s >= 0]

        TN_orig = np.zeros((len(preserve_ist_list), len(orig_classes_with_aux)))
        for pi, ist in enumerate(preserve_ist_list):
            if ist < TN.shape[0]:
                for oi, orig_class in enumerate(orig_classes_with_aux):
                    if orig_class < TN.shape[1]:
                        TN_orig[pi, oi] = TN[ist, orig_class]

        # MATLAB lines 306-319: Merge auxiliary classes
        for r_idx in range(len(fjclassmap)):
            s = int(fjclassmap[r_idx])
            if s >= 0 and s < orig_nclasses:
                aux_col = orig_nclasses + r_idx
                if aux_col < QN.shape[1]:
                    QN[:, s] = QN[:, s] + QN[:, aux_col]
                    UN[:, s] = UN[:, s] + UN[:, aux_col]
                    TN[:, s] = TN[:, s] + TN[:, aux_col]
                    # RN = QN / TN (MATLAB line 316)
                    with np.errstate(divide='ignore', invalid='ignore'):
                        RN[:, s] = np.where(TN[:, s] > 0, QN[:, s] / TN[:, s], 0)

        # MATLAB line 321: Restore TN at Join and Source
        for pi, ist in enumerate(preserve_ist_list):
            if ist < TN.shape[0]:
                for oi, orig_class in enumerate(orig_classes_with_aux):
                    if orig_class < TN.shape[1]:
                        TN[ist, orig_class] = TN_orig[pi, oi]

        # Zero out UN at Join stations (Join has Immediate service, so utilization must be 0)
        for j_idx in join_indices:
            if j_idx < len(nonfj_sn.nodeToStation):
                join_ist = int(nonfj_sn.nodeToStation[j_idx])
                if join_ist < UN.shape[0]:
                    UN[join_ist, :orig_nclasses] = 0

        # MATLAB lines 323-328: Delete auxiliary class columns
        QN_final = QN[:orig_nstations, :orig_nclasses].copy()
        UN_final = UN[:orig_nstations, :orig_nclasses].copy()
        RN_final = RN[:orig_nstations, :orig_nclasses].copy()
        TN_final = TN[:orig_nstations, :orig_nclasses].copy()
        CN_final = CN[:orig_nclasses].copy() if len(CN) >= orig_nclasses else np.sum(RN_final, axis=0)
        XN_final = XN[:orig_nclasses].copy() if len(XN) >= orig_nclasses else np.zeros(orig_nclasses)

        # Compute arrival rates using original sn (MATLAB line 336)
        from ..api.sn.getters import sn_get_arvr_from_tput
        AN_final = sn_get_arvr_from_tput(orig_sn, TN_final, None)

        # Compute residence times (MATLAB line 337)
        WN_final = sn_get_residt_from_respt(orig_sn, RN_final, None)

        self._result = self._fj_publish({
            'QN': QN_final,
            'UN': UN_final,
            'RN': RN_final,
            'TN': TN_final,
            'AN': AN_final,
            'XN': XN_final,
            'WN': WN_final,
            'CN': np.sum(RN_final, axis=0),
            'runtime': 0.0,
            'method': 'mmt_open',
            'iter': total_iter
        })

        return self._result

    @staticmethod
    def _combine_routing_matrix(P):
        """Combine cell array routing matrix into single matrix."""
        if P is None:
            return None
        if isinstance(P, np.ndarray):
            return P
        # P is list of lists
        nclasses = len(P)
        if nclasses == 0:
            return None
        nnodes = P[0][0].shape[0] if P[0][0] is not None else 0
        combined = np.zeros((nclasses * nnodes, nclasses * nnodes))
        for r in range(nclasses):
            for s in range(nclasses):
                if r < len(P) and s < len(P[r]) and P[r][s] is not None:
                    combined[r*nnodes:(r+1)*nnodes, s*nnodes:(s+1)*nnodes] = P[r][s]
        return combined

    def _run_ht_transformation(self):
        """
        Run full Heidelberger-Trivedi transformation for closed fork-join networks.

        Implements the full H-T method from MATLAB:
        1. Transform model ONCE using ModelAdapter.ht() to create auxiliary classes
        2. Iterate:
           a. Run standard MVA on the transformed model
           b. Check convergence on QN
           c. Update sync delays at Join and AuxDelay nodes
        3. Merge auxiliary class results back into original classes

        Returns:
            Result dictionary or None if transformation fails
        """
        from ..io.model_adapter import ModelAdapter
        from ..api.sn.network_struct import NodeType
        from ..distributions.continuous import Exp
        from itertools import combinations

        sn = self._sn
        if sn is None or sn.fj is None or not np.any(sn.fj):
            return None

        # Step 1: Transform model using ModelAdapter.ht() (ONCE at the start)
        try:
            ht_result = ModelAdapter.ht(self.model)
        except Exception as e:
            import warnings
            warnings.warn(f"H-T transformation failed: {e}. Falling back to approximate method.")
            return self._run_approximate_ht()

        nonfjmodel = ht_result.nonfjmodel
        fjclassmap = ht_result.fjclassmap
        fjforkmap = ht_result.fjforkmap
        fj_auxiliary_delays = ht_result.fj_auxiliary_delays

        import os
        debug_ht = os.environ.get('DEBUG_HT', '0') == '1'

        if debug_ht:
            print("\n=== H-T Transformation Results ===")
            print(f"Original model: {len(self.model._classes)} classes, {len(self.model._nodes)} nodes")
            print(f"Transformed model: {len(nonfjmodel._classes)} classes, {len(nonfjmodel._nodes)} nodes")
            print(f"fjclassmap: {fjclassmap}")
            print(f"fjforkmap: {fjforkmap}")
            print(f"fj_auxiliary_delays: {fj_auxiliary_delays}")

        # Check if transformation actually created auxiliary classes
        if len(nonfjmodel._classes) == len(self.model._classes):
            if debug_ht:
                print("No auxiliary classes created, falling back to approximate H-T")
            return self._run_approximate_ht()

        orig_sn = sn
        orig_nclasses = len(self.model._classes)

        # Find fork indices in original model
        fork_indices = []
        for i, nt in enumerate(orig_sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK.value if hasattr(NodeType.FORK, 'value') else NodeType.FORK:
                fork_indices.append(i)

        # Iteration parameters
        max_iter = self.options.max_iter if hasattr(self.options, 'max_iter') else 1000
        coarse_tol = 1e-3

        # Initialize for convergence check
        QN_prev = None
        converged = False
        iter_count = 0

        if debug_ht:
            print(f"\nFork indices: {fork_indices}")
            print(f"Max iterations: {max_iter}, Tolerance: {coarse_tol}")

        while not converged and iter_count < max_iter:
            iter_count += 1

            # Refresh and run MVA on transformed model
            nonfjmodel.refresh_struct()
            nonfjstruct = nonfjmodel._sn

            if debug_ht and iter_count == 1:
                print("\n=== Transformed Model Structure ===")
                print(f"Stations: {nonfjstruct.nstations}, Classes: {nonfjstruct.nclasses}, Chains: {nonfjstruct.nchains}")
                print(f"stationToNode: {nonfjstruct.stationToNode}")
                print(f"nodeToStation: {nonfjstruct.nodeToStation}")
                print(f"refstat: {nonfjstruct.refstat}")
                print(f"njobs (populations): {nonfjstruct.njobs}")
                print(f"chains: {nonfjstruct.chains}")
                print("Service rates (rates matrix):")
                print(nonfjstruct.rates)
                print("Visits per chain:")
                for c in range(len(nonfjstruct.visits)):
                    print(f"  Chain {c}: {nonfjstruct.visits[c]}")

            nonfj_solver = self._fj_inner_solver(nonfjmodel, method=getattr(self, 'method', None))
            nonfj_solver._skip_fork_join = True

            try:
                nonfj_solver.runAnalyzer()
            except Exception as e:
                import warnings
                warnings.warn(f"MVA on transformed model failed: {e}. Falling back to approximate method.")
                return self._run_approximate_ht()

            result = self._fj_result_dict(nonfj_solver)
            if result is None:
                return self._run_approximate_ht()
            QN = result.get('QN', np.zeros((nonfj_solver.nstations, nonfj_solver.nclasses)))
            RN = result.get('RN', np.zeros((nonfj_solver.nstations, nonfj_solver.nclasses)))
            TN = result.get('TN', np.zeros((nonfj_solver.nstations, nonfj_solver.nclasses)))
            UN = result.get('UN', np.zeros((nonfj_solver.nstations, nonfj_solver.nclasses)))
            XN = result.get('XN', np.zeros(nonfj_solver.nclasses))

            if debug_ht:
                print(f"\n=== Iteration {iter_count} ===")
                qn_per_class = np.nansum(QN, axis=0)
                print(f"QN sum per class: {qn_per_class}")
                print(f"XN: {XN}")
                print(f"RN matrix:\n{RN}")

            # Check convergence (MATLAB: max(abs(1-QN_1./QN)) < CoarseTol && forkIter > 2)
            max_diff = float('inf')
            if QN_prev is not None and iter_count > 2:
                # Use full QN matrix for convergence check (like MATLAB)
                # Only consider entries where QN is non-negligible (to avoid 0/0 issues)
                from ..constants import GlobalConstants
                # Create mask for meaningful entries (both current and previous are > Zero)
                meaningful = (np.abs(QN) > GlobalConstants.Zero) & (np.abs(QN_prev) > GlobalConstants.Zero)
                if np.any(meaningful):
                    ratio = np.divide(QN_prev, QN, out=np.ones_like(QN), where=meaningful)
                    diff = np.abs(1 - ratio)
                    diff = np.where(meaningful, diff, 0.0)  # Ignore non-meaningful entries
                    max_diff = np.max(diff)
                else:
                    max_diff = 0.0  # All zeros, consider converged
                if debug_ht:
                    print(f"  Convergence check: max_diff={max_diff:.6e}, tol={coarse_tol}")
                if max_diff < coarse_tol:
                    converged = True
                    if debug_ht:
                        print(f"  CONVERGED at iteration {iter_count}")

            QN_prev = QN.copy()

            if converged:
                break

            # Update sync delays (MATLAB lines 238-272)
            from ..constants import GlobalConstants

            if debug_ht:
                print(f"\n--- Sync Delay Update (iter {iter_count}) ---")

            for f in fork_indices:
                join_idx_arr = np.where(orig_sn.fj[f, :] > 0)[0]
                if len(join_idx_arr) == 0:
                    continue
                join_idx = join_idx_arr[0]

                aux_delay_idx = fj_auxiliary_delays.get(join_idx, None)
                if aux_delay_idx is None:
                    if debug_ht:
                        print(f"  Fork {f}: No aux_delay for join_idx={join_idx}")
                    continue

                # Get station indices in transformed model
                join_station = nonfjstruct.nodeToStation[join_idx] if join_idx < len(nonfjstruct.nodeToStation) else -1
                aux_delay_station = nonfjstruct.nodeToStation[aux_delay_idx] if aux_delay_idx < len(nonfjstruct.nodeToStation) else -1

                if debug_ht:
                    print(f"  Fork {f}: join_idx={join_idx}, join_station={join_station}, aux_delay_idx={aux_delay_idx}, aux_delay_station={aux_delay_station}")

                if join_station < 0 or aux_delay_station < 0:
                    continue

                # Get tasksPerLink (fanout) from fork node
                # MATLAB: sn.nodeparam{f}.fanOut (tasksPerLink, typically 1)
                fork_node = self.model._nodes[f]
                fanout = 1  # Default: 1 task per forked branch
                if hasattr(fork_node, 'get_tasks_per_link'):
                    tpl = fork_node.get_tasks_per_link()
                    if tpl is not None and not np.any(np.isnan(tpl)):
                        fanout = int(np.max(tpl)) if hasattr(tpl, '__len__') else int(tpl)
                        if fanout < 1:
                            fanout = 1

                if debug_ht:
                    print(f"    fanout (tasksPerLink): {fanout}")

                # For each original class
                for r in range(orig_nclasses):
                    # Check if this class visits the fork (MATLAB line 243: if sn.nodevisits{c}(f,r) == 0 continue)
                    # We need to check nodevisits to skip classes that don't use this fork
                    skip_class = False
                    if hasattr(orig_sn, 'nodevisits') and orig_sn.nodevisits is not None:
                        for c in range(orig_sn.nchains):
                            if hasattr(orig_sn, 'chains') and orig_sn.chains is not None:
                                chains_arr = np.asarray(orig_sn.chains)
                                if chains_arr.ndim == 2 and c < chains_arr.shape[0] and r < chains_arr.shape[1]:
                                    if chains_arr[c, r] > 0:  # Class r is in chain c
                                        nodevisits = orig_sn.nodevisits[c]
                                        if nodevisits is not None and f < nodevisits.shape[0] and r < nodevisits.shape[1]:
                                            if nodevisits[f, r] == 0:
                                                skip_class = True
                                                break

                    if skip_class:
                        if debug_ht:
                            print(f"    Class {r}: Skipping (nodevisits[{f},{r}] == 0)")
                        continue

                    # Get auxiliary classes for this original class and fork
                    aux_class_indices = np.where((fjclassmap == r) & (fjforkmap == f))[0]
                    if len(aux_class_indices) == 0:
                        if debug_ht:
                            print(f"    Class {r}: No auxiliary classes")
                        continue

                    if debug_ht:
                        print(f"    Class {r}: aux_class_indices={aux_class_indices}")

                    # MATLAB: ri = RN(:, find(fjclassmap == r));
                    # MATLAB: ri = sum(ri, 1, "omitnan") - RN(aux_delay, aux_classes) - RN(join, aux_classes)
                    ri = []
                    for aux_idx in aux_class_indices:
                        # Sum RN across all stations for this aux class
                        total_rn = np.nansum(RN[:, aux_idx])
                        # Subtract join and aux delay
                        rn_join = RN[join_station, aux_idx] if not np.isnan(RN[join_station, aux_idx]) else 0
                        rn_auxdelay = RN[aux_delay_station, aux_idx] if not np.isnan(RN[aux_delay_station, aux_idx]) else 0
                        branch_rt = total_rn - rn_join - rn_auxdelay
                        if debug_ht:
                            print(f"      aux_idx={aux_idx}: sum(RN[:,aux])={total_rn:.6f}, RN[join,aux]={rn_join:.6f}, RN[auxdelay,aux]={rn_auxdelay:.6f} => ri={branch_rt:.6f}")
                        # Handle NaN/Inf
                        if np.isnan(branch_rt) or np.isinf(branch_rt):
                            branch_rt = 0
                        ri.append(max(branch_rt, 1e-10))  # Avoid division by zero

                    ri_arr = np.array(ri)
                    num_branches = len(ri_arr)

                    if debug_ht:
                        print(f"    ri_arr={ri_arr}, num_branches={num_branches}")

                    if num_branches < 2:
                        if debug_ht:
                            print(f"    Skipping: num_branches < 2")
                        continue

                    # MATLAB line 252: parallel_branches = length(self.model.nodes{f}.output.outputStrategy{r}{3});
                    # Get the actual number of parallel branches from fork output strategy
                    parallel_branches = num_branches
                    if hasattr(fork_node, 'output') and hasattr(fork_node.output, 'output_strategy'):
                        strategy = fork_node.output.output_strategy
                        if isinstance(strategy, dict) and r in strategy:
                            strat_r = strategy[r]
                            if isinstance(strat_r, (list, tuple)) and len(strat_r) >= 3:
                                parallel_branches = len(strat_r[2]) if hasattr(strat_r[2], '__len__') else num_branches
                    if debug_ht:
                        print(f"    parallel_branches={parallel_branches}")

                    # Compute E[max] using H-T inclusion-exclusion formula
                    lambdai = 1.0 / ri_arr
                    d0 = 0.0
                    for pow_val in range(parallel_branches):
                        combos = list(combinations(lambdai, pow_val + 1))
                        current_sum = sum(1.0 / sum(combo) for combo in combos)
                        d0 += ((-1) ** pow_val) * current_sum

                    if debug_ht:
                        print(f"    lambdai={lambdai}, d0={d0:.6f}")

                    # Individual sync delays: di = d0 * fanout - ri
                    di = d0 * fanout - ri_arr

                    if debug_ht:
                        print(f"    di = d0*fanout - ri = {d0:.6f}*{fanout} - {ri_arr} = {di}")

                    # r0 = sum of response times for all inchain classes, minus join
                    # MATLAB: r0 = sum(RN(:, inchain), 2); r0 = sum(r0, 1, "omitnan") - RN(joinIdx, r)
                    # Find inchain classes (all classes in the same chain as r)
                    inchain = [r]  # Start with original class
                    if hasattr(nonfjstruct, 'inchain') and nonfjstruct.inchain is not None:
                        chain_id = get_chain_for_class(nonfjstruct.chains, r)
                        if chain_id >= 0 and chain_id in nonfjstruct.inchain:
                            inchain = list(nonfjstruct.inchain[chain_id])
                    elif hasattr(nonfjstruct, 'chains') and nonfjstruct.chains is not None:
                        chain_id = get_chain_for_class(nonfjstruct.chains, r)
                        if chain_id >= 0:
                            # Fallback: find classes by searching 2D chains matrix
                            chains_arr = np.asarray(nonfjstruct.chains)
                            if chains_arr.ndim == 2 and chain_id < chains_arr.shape[0]:
                                inchain = list(np.where(chains_arr[chain_id, :] > 0)[0])

                    if debug_ht:
                        print(f"    inchain={inchain}")

                    # Sum RN across all inchain classes
                    # MATLAB: r0 = sum(RN(:, inchain), 2) gives a column vector, then sum(..., 1) sums it
                    r0_matrix = RN[:, inchain]  # M x len(inchain)
                    r0_matrix = np.where(np.isnan(r0_matrix) | np.isinf(r0_matrix), 0, r0_matrix)
                    r0 = np.sum(r0_matrix) - RN[join_station, r]
                    if np.isnan(r0) or np.isinf(r0):
                        r0 = 0
                    r0 = max(r0, GlobalConstants.FineTol)

                    if debug_ht:
                        print(f"    r0 = sum(RN[:, inchain]) - RN[join, r] = {np.sum(r0_matrix):.6f} - {RN[join_station, r]:.6f} = {r0:.6f}")

                    # Update service times
                    # MATLAB line 262: Join service for original class r: d0 * fanout
                    sync_delay_orig = d0 * fanout
                    if debug_ht:
                        print(f"    Setting join[{join_idx}].service(class {r}) = d0*fanout = {sync_delay_orig:.6f}")
                    if sync_delay_orig > 1e-10:
                        nonfjmodel._nodes[join_idx].set_service(nonfjmodel._classes[r], Exp.fit_mean(sync_delay_orig))

                    # For each auxiliary class (MATLAB lines 264-268)
                    for idx, aux_idx in enumerate(aux_class_indices):
                        # Join service for aux class: di[idx]
                        aux_sync_delay = max(di[idx], 1e-10)
                        if debug_ht:
                            print(f"    Setting join[{join_idx}].service(aux_class {aux_idx}) = di[{idx}] = {aux_sync_delay:.6f}")
                            print(f"    Setting aux_delay[{aux_delay_idx}].service(aux_class {aux_idx}) = r0 = {r0:.6f}")
                        nonfjmodel._nodes[join_idx].set_service(nonfjmodel._classes[aux_idx], Exp.fit_mean(aux_sync_delay))
                        # Aux delay service for aux class: r0
                        nonfjmodel._nodes[aux_delay_idx].set_service(nonfjmodel._classes[aux_idx], Exp.fit_mean(r0))

        # Step 3: Merge results (MATLAB lines 274-302)
        if debug_ht:
            print(f"\n=== Step 3: Merge Results ===")
            print(f"Final QN before merge:\n{QN}")
            print(f"Final RN before merge:\n{RN}")
            print(f"Final TN before merge:\n{TN}")
            print(f"Final XN before merge: {XN}")

        # Find join and aux delay station indices
        join_stations = []
        aux_delay_stations = []
        for f in fork_indices:
            join_idx_arr = np.where(orig_sn.fj[f, :] > 0)[0]
            if len(join_idx_arr) > 0:
                join_idx = join_idx_arr[0]
                if join_idx < len(nonfjstruct.nodeToStation):
                    join_stations.append(nonfjstruct.nodeToStation[join_idx])
                aux_delay_idx = fj_auxiliary_delays.get(join_idx, None)
                if aux_delay_idx is not None and aux_delay_idx < len(nonfjstruct.nodeToStation):
                    aux_delay_stations.append(nonfjstruct.nodeToStation[aux_delay_idx])

        # Get original class indices that have auxiliary classes
        orig_classes_with_aux = sorted(set(fjclassmap[fjclassmap >= 0]))

        if debug_ht:
            print(f"join_stations: {join_stations}")
            print(f"aux_delay_stations: {aux_delay_stations}")
            print(f"orig_classes_with_aux: {orig_classes_with_aux}")
            print(f"fjclassmap: {fjclassmap}")
            print(f"stations_to_keep will be: {[i for i in range(QN.shape[0]) if i not in aux_delay_stations]}")

        # Save original class throughputs at join stations before zeroing
        TN_orig_join = {}
        for js in join_stations:
            for oc in orig_classes_with_aux:
                if js < TN.shape[0] and oc < TN.shape[1]:
                    TN_orig_join[(js, oc)] = TN[js, oc]

        if debug_ht:
            print(f"TN_orig_join: {TN_orig_join}")

        # Zero out original class metrics at join stations
        for js in join_stations:
            for oc in orig_classes_with_aux:
                if js < QN.shape[0] and oc < QN.shape[1]:
                    QN[js, oc] = 0
                    RN[js, oc] = 0
                    TN[js, oc] = 0
                    UN[js, oc] = 0

        # Remove auxiliary delay station rows (create mask for stations to keep)
        stations_to_keep = [i for i in range(QN.shape[0]) if i not in aux_delay_stations]

        QN_merged = QN[stations_to_keep, :]
        UN_merged = UN[stations_to_keep, :]
        RN_merged = RN[stations_to_keep, :]
        TN_merged = TN[stations_to_keep, :]

        if debug_ht:
            print(f"After removing aux delay stations:")
            print(f"QN_merged shape: {QN_merged.shape}")
            print(f"QN_merged:\n{QN_merged}")
            print(f"RN_merged:\n{RN_merged}")
            print(f"TN_merged:\n{TN_merged}")

        # Merge auxiliary classes into original classes
        # fjclassmap[class_idx] = -1 for original classes, = orig_class_idx for auxiliary classes
        # Auxiliary classes are stored at their actual indices in the fjclassmap
        for aux_class_idx in range(len(fjclassmap)):
            orig_idx = fjclassmap[aux_class_idx]
            # Only process auxiliary classes (those that map to an original class)
            if orig_idx >= 0 and aux_class_idx < QN_merged.shape[1] and orig_idx < QN_merged.shape[1]:
                if debug_ht:
                    print(f"Merging aux_class {aux_class_idx} into orig_class {orig_idx}:")
                    print(f"  QN_merged[:,{orig_idx}] before: {QN_merged[:, orig_idx]}")
                    print(f"  QN_merged[:,{aux_class_idx}]: {QN_merged[:, aux_class_idx]}")
                QN_merged[:, orig_idx] += QN_merged[:, aux_class_idx]
                UN_merged[:, orig_idx] += UN_merged[:, aux_class_idx]
                # Add all throughputs of the auxiliary classes to facilitate the computation of the response times
                TN_merged[:, orig_idx] += TN_merged[:, aux_class_idx]
                if debug_ht:
                    print(f"  QN_merged[:,{orig_idx}] after: {QN_merged[:, orig_idx]}")
                # Recompute RN = QN / TN after merging (matching JAR behavior)
                for i in range(RN_merged.shape[0]):
                    if TN_merged[i, orig_idx] > 0:
                        RN_merged[i, orig_idx] = QN_merged[i, orig_idx] / TN_merged[i, orig_idx]
                if debug_ht:
                    print(f"  RN_merged[:,{orig_idx}] after: {RN_merged[:, orig_idx]}")

        if debug_ht:
            print(f"After merging auxiliary classes:")
            print(f"QN_merged:\n{QN_merged}")
            print(f"RN_merged:\n{RN_merged}")
            print(f"TN_merged:\n{TN_merged}")

        # Restore original class throughputs at join stations
        for (old_js, oc), tput in TN_orig_join.items():
            if old_js in stations_to_keep:
                new_js_idx = stations_to_keep.index(old_js)
                if new_js_idx < TN_merged.shape[0] and oc < TN_merged.shape[1]:
                    TN_merged[new_js_idx, oc] = tput

        # Zero out UN at Join stations for all original classes
        # Join has Immediate service in the original model, so utilization must be 0
        for js in join_stations:
            if js in stations_to_keep:
                new_js_idx = stations_to_keep.index(js)
                if new_js_idx < UN_merged.shape[0]:
                    UN_merged[new_js_idx, :orig_nclasses] = 0

        # Extract only original classes
        QN_final = QN_merged[:, :orig_nclasses]
        UN_final = UN_merged[:, :orig_nclasses]
        RN_final = RN_merged[:, :orig_nclasses]
        TN_final = TN_merged[:, :orig_nclasses]
        XN_final = XN[:orig_nclasses].copy()

        if debug_ht:
            print(f"\n=== Final Results (orig_nclasses={orig_nclasses}) ===")
            print(f"QN_final:\n{QN_final}")
            print(f"RN_final:\n{RN_final}")
            print(f"TN_final:\n{TN_final}")
            print(f"XN_final: {XN_final}")
            print(f"CN (sum RN): {np.sum(RN_final, axis=0)}")

        # Compute arrival rates using sn_get_arvr_from_tput (MATLAB runAnalyzer.m line 341)
        from ..api.sn.getters import sn_get_arvr_from_tput
        AN_final = sn_get_arvr_from_tput(self._sn, TN_final, None)

        # Compute residence times from response times (WN = RN * V)
        from ..api.sn.transforms import sn_get_residt_from_respt
        WN_final = sn_get_residt_from_respt(self._sn, RN_final, None)

        # Store results
        self._result = self._fj_publish({
            'QN': QN_final,
            'UN': UN_final,
            'RN': RN_final,
            'TN': TN_final,
            'AN': AN_final,
            'XN': XN_final,
            'WN': WN_final,
            'CN': np.sum(RN_final, axis=0),
            'runtime': 0.0,
            'method': 'ht',
            'iter': iter_count
        })

        return self._result

    def _run_approximate_ht(self):
        """
        Run Heidelberger-Trivedi MVA for closed fork-join networks.

        Implements the H-T approximation which:
        1. Analyzes each parallel branch independently with full population N
        2. Computes E[max(R_1, ..., R_K)] for synchronization delay
        3. Uses cycle time C = Z + E[max] to compute throughput X = N/C
        4. Iterates until queue lengths converge (matching MATLAB's approach)

        The H-T formula uses exponential approximation for the maximum:
        E[max] = sum_{k=1}^{K} (-1)^{k+1} * sum_{combos of k} 1/sum(rates)

        Note: Results typically differ from exact solutions by ~4% due to the
        exponential assumption in the sync delay formula.

        Returns:
            Result dictionary or None if analysis fails
        """
        from ..api.sn.network_struct import NodeType

        sn = self._sn
        if sn is None or sn.fj is None or not np.any(sn.fj):
            return None

        # Find fork and join indices
        fork_indices = []
        join_indices = []
        for i, nt in enumerate(sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK.value if hasattr(NodeType.FORK, 'value') else NodeType.FORK:
                fork_indices.append(i)
            elif nt_val == NodeType.JOIN.value if hasattr(NodeType.JOIN, 'value') else NodeType.JOIN:
                join_indices.append(i)

        if len(fork_indices) == 0:
            return None

        # Get queueing demands (service times at queues) - also returns queue indices
        L, queue_indices = self._get_queueing_demands()

        # Get delay station indices
        delay_indices = self._get_delay_stations()

        # Get think times (delay station service times)
        Z_base = self._get_think_times()

        # Handle pure delay network (no queueing stations)
        # This happens when all stations are delays or INF (infinite servers)
        if L.shape[0] == 0:
            R = self.nclasses
            QN = np.zeros((self.nstations, R))
            UN = np.zeros((self.nstations, R))
            RN = np.zeros((self.nstations, R))
            TN = np.zeros((self.nstations, R))
            AN = np.zeros((self.nstations, R))
            XN = np.zeros(R)

            # For pure delay network:
            # X[r] = N[r] / Z[r] (throughput = population / think time)
            # R[i,r] = D[i,r] (response time = service time, no queueing)
            # Q[i,r] = X[r] * D[i,r] (using Little's law)
            N = self.njobs.copy()
            for r in range(R):
                if np.isfinite(N[r]) and N[r] > 0 and Z_base[r] > 0:
                    XN[r] = N[r] / Z_base[r]
                    for i in range(self.nstations):
                        if self.demands[i, r] > 0:
                            RN[i, r] = self.demands[i, r]  # Response = service time
                            TN[i, r] = XN[r]  # Throughput at each station
                            QN[i, r] = XN[r] * self.demands[i, r]  # Queue length (Little's law)
                            UN[i, r] = QN[i, r]  # Utilization = queue length for delay stations

            return QN, UN, RN, TN, AN, XN

        # Get population
        N = self.njobs.copy()

        # Initialize result arrays
        R = self.nclasses
        QN = np.zeros((self.nstations, R))
        UN = np.zeros((self.nstations, R))
        RN = np.zeros((self.nstations, R))
        TN = np.zeros((self.nstations, R))
        AN = np.zeros((self.nstations, R))
        XN = np.zeros(R)

        # Identify parallel branches for each fork-join pair
        # Use branch structure to properly handle asymmetric branches (serial stations within a branch)
        fj_branches = {}
        for f in fork_indices:
            join_idx_arr = np.where(sn.fj[f, :])[0]
            if len(join_idx_arr) > 0:
                join_idx = join_idx_arr[0]
                join_station = sn.nodeToStation[join_idx] if join_idx < len(sn.nodeToStation) else -1
                # Get branch structure (list of branches, each branch is list of stations)
                branches = self._find_branch_structure(f, join_idx)
                # Also keep flat list for backward compatibility
                parallel_stations = self._find_parallel_branches(f, join_idx)
                fj_branches[f] = {
                    'join_node': join_idx,
                    'join_station': join_station,
                    'parallel_stations': parallel_stations,
                    'branches': branches  # Each branch is a list of serial stations
                }

        # Initial response times (just service times)
        R_branch = {}  # Response times at parallel branch queues
        for f, fj_info in fj_branches.items():
            R_branch[f] = {}
            for r in range(R):
                if np.isfinite(N[r]) and N[r] > 0:
                    branch_ri = []
                    for station_idx in fj_info['parallel_stations']:
                        if station_idx in queue_indices:
                            q_idx = queue_indices.index(station_idx)
                            branch_ri.append(L[q_idx, r])  # Initial: service time
                        else:
                            branch_ri.append(self.demands[station_idx, r])
                    R_branch[f][r] = np.array(branch_ri) if branch_ri else np.array([1.0])

        # H-T MVA for closed fork-join with outer convergence loop
        # This matches MATLAB's approach which iterates until QN converges
        # Key H-T insight: Each branch is analyzed INDEPENDENTLY with full population N
        max_outer_iter = self.options.max_iter if hasattr(self.options, 'max_iter') else 1000
        coarse_tol = 1e-3  # Convergence tolerance (MATLAB's CoarseTol)

        # Initialize sync delays to 0
        sync_delays = {r: 0.0 for r in range(R)}
        QN_prev = np.ones((self.nstations, R)) * 1e10  # Large initial value for convergence check

        outer_iter = 0
        converged = False

        while not converged and outer_iter < max_outer_iter:
            outer_iter += 1

            # For each class, run H-T analysis
            for r in range(R):
                if not (np.isfinite(N[r]) and N[r] > 0):
                    continue

                n_jobs = int(N[r])

                for f, fj_info in fj_branches.items():
                    branches = fj_info['branches']
                    num_branches = len(branches)

                    if num_branches == 0:
                        continue

                    # H-T approach: analyze each branch INDEPENDENTLY with full population N
                    # Each branch is treated as a separate closed network
                    branch_R = np.zeros(num_branches)  # Branch response times
                    branch_Q = {}  # Per-station queue lengths

                    for b, branch in enumerate(branches):
                        # Get demands for stations in this branch
                        branch_demands = []
                        branch_stations = []
                        for station_idx in branch:
                            if self.rates[station_idx, r] > 0:
                                D_s = 1.0 / self.rates[station_idx, r]
                                branch_demands.append(D_s)
                                branch_stations.append(station_idx)

                        if len(branch_stations) == 0:
                            branch_R[b] = 1e-10
                            continue

                        # Run MVA for this branch alone with population N
                        # Initial queue lengths are 0
                        station_Q = {s: 0.0 for s in branch_stations}
                        station_R = {s: branch_demands[i] for i, s in enumerate(branch_stations)}

                        for k in range(1, n_jobs + 1):
                            # Response time at each station: R_i = D_i * (1 + Q_i)
                            for i, s in enumerate(branch_stations):
                                D_s = branch_demands[i]
                                station_R[s] = D_s * (1 + station_Q[s])

                            # Total branch response time
                            R_branch_k = sum(station_R[s] for s in branch_stations)

                            # Cycle time for this branch = Z + R_branch + sync_delay
                            # (sync_delay accounts for waiting for other branches from previous iter)
                            C_k = Z_base[r] + R_branch_k + sync_delays[r]

                            # Throughput for k jobs
                            X_k = k / C_k

                            # Update queue lengths
                            for s in branch_stations:
                                station_Q[s] = X_k * station_R[s]

                        # Store final branch response time and queue lengths
                        branch_R[b] = sum(station_R[s] for s in branch_stations)
                        for s in branch_stations:
                            branch_Q[s] = station_Q[s]

                    # Compute E[max] of branch response times
                    if num_branches >= 2:
                        d0 = self._compute_sync_delay(branch_R)
                    else:
                        d0 = branch_R[0]

                    # Final cycle time and throughput
                    C_final = Z_base[r] + d0
                    X_final = n_jobs / C_final

                    # Store final results
                    XN[r] = X_final
                    R_branch[f][r] = branch_R

                    # Store individual station results
                    # In H-T, each auxiliary class visits only one branch with full population N
                    # The queue length at each station is from its auxiliary class
                    # Response time RN = QN / TN (Little's Law)
                    for b, branch in enumerate(branches):
                        for station_idx in branch:
                            if station_idx in branch_Q:
                                D_s = 1.0 / self.rates[station_idx, r] if self.rates[station_idx, r] > 0 else 0
                                # Queue length is from the auxiliary class (not divided by branches)
                                QN[station_idx, r] = branch_Q[station_idx]
                                # Utilization = X * D for the actual system throughput
                                UN[station_idx, r] = X_final * D_s
                                TN[station_idx, r] = X_final
                                # Response time = QN / TN (Little's Law)
                                RN[station_idx, r] = QN[station_idx, r] / TN[station_idx, r] if TN[station_idx, r] > 0 else 0
                                # Arrival rate at parallel stations = throughput / num_branches
                                AN[station_idx, r] = X_final / num_branches

                    # Update sync delay for next iteration
                    sync_delays[r] = max(0, d0 - np.mean(branch_R))

            # Set delay station metrics - use individual station demands, not total Z
            for station_idx in delay_indices:
                for r in range(R):
                    D_i = self.demands[station_idx, r]
                    if D_i > 0:
                        RN[station_idx, r] = D_i
                        QN[station_idx, r] = XN[r] * D_i
                        UN[station_idx, r] = XN[r] * D_i  # Utilization at delay
                        TN[station_idx, r] = XN[r]
                        AN[station_idx, r] = XN[r]

            # Set join station metrics
            for f, fj_info in fj_branches.items():
                join_station = fj_info['join_station']
                if join_station >= 0 and join_station < self.nstations:
                    for r in range(R):
                        if np.isfinite(N[r]) and N[r] > 0:
                            ri = R_branch[f].get(r, np.array([1.0]))
                            if len(ri) >= 2:
                                d0 = self._compute_sync_delay(ri)
                                sync_d = max(0, d0 - np.mean(ri))
                                num_branches = len(ri)
                            else:
                                sync_d = 0
                                num_branches = 1

                            RN[join_station, r] = sync_d
                            QN[join_station, r] = XN[r] * sync_d * num_branches
                            TN[join_station, r] = XN[r]
                            # Arrival rate at join = throughput * num_branches (jobs from all branches merge)
                            AN[join_station, r] = XN[r] * num_branches
                            UN[join_station, r] = 0  # Join has no service

            # Check convergence: max(abs(1 - QN_prev / QN)) < tol
            # Avoid division by zero
            QN_safe = np.where(QN > 0, QN, 1e-10)
            QN_prev_safe = np.where(QN_prev > 0, QN_prev, 1e-10)
            rel_change = np.max(np.abs(1 - QN_prev_safe / QN_safe))

            if outer_iter >= 2 and rel_change < coarse_tol:
                converged = True

            QN_prev = QN.copy()

        # Compute final cycle times
        CN = np.zeros(R)
        for r in range(R):
            if np.isfinite(N[r]) and N[r] > 0 and XN[r] > 0:
                CN[r] = N[r] / XN[r]
            else:
                CN[r] = 0

        # Compute residence times from response times (WN = RN * V)
        from ..api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self._sn, RN, None)

        # Store results
        self._result = self._fj_publish({
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
            'WN': WN,
            'CN': CN,
            'method': 'approximate_ht',
            'iter': outer_iter
        })

        return self._result

    def _compute_fork_join_sync_delays(self, QN, UN, RN, TN, AN, XN):
        """
        Compute and add synchronization delays at Join nodes.

        For each fork-join pair, computes the expected synchronization delay
        using the Heidelberger-Trivedi formula.

        Note: This is a lightweight fallback used only when
        ``_run_fork_join_analysis`` (the primary path, which ports MATLAB's full
        MMT/H-T model transformation via ModelAdapter.mmt/ht) returns None. On the
        primary path native MVA matches MATLAB on closed fork-join networks; this
        post-processing sync-delay estimate is a coarser approximation retained
        only as a safety net.

        Args:
            QN, UN, RN, TN, AN, XN: Performance metric arrays (modified in place)
        """
        from ..api.sn.network_struct import NodeType
        from itertools import combinations

        sn = self._sn
        if sn.fj is None or not np.any(sn.fj):
            return

        # Find fork indices
        fork_indices = []
        for i, nt in enumerate(sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK:
                fork_indices.append(i)

        for f in fork_indices:
            # Find join associated with this fork
            join_idx_arr = np.where(sn.fj[f, :])[0]
            if len(join_idx_arr) == 0:
                continue
            join_idx = join_idx_arr[0]

            # Get join station index
            if join_idx >= len(sn.nodeToStation):
                continue
            join_station = sn.nodeToStation[join_idx]
            if join_station < 0 or join_station >= self.nstations:
                continue

            # Find the branch structure (which stations are on each branch)
            branches = self._find_branch_structure(f, join_idx)
            parallel_stations = self._find_parallel_branches(f, join_idx)

            if len(branches) < 2:
                continue

            # Get branch response times for each class
            # Sum response times for serial stations within each branch
            for r in range(self.nclasses):
                ri = []
                for branch in branches:
                    branch_resp_time = 0.0
                    for station_idx in branch:
                        if RN[station_idx, r] > 0:
                            branch_resp_time += RN[station_idx, r]
                    if branch_resp_time > 0:
                        ri.append(branch_resp_time)

                if len(ri) >= 2:
                    ri_arr = np.array(ri)
                    num_branches = len(ri_arr)

                    # Compute expected max using H-T formula with response times
                    d0 = self._compute_sync_delay(ri_arr)
                    # Sync delay at join = E[max] - mean (waiting for slowest)
                    sync_delay = d0 - np.mean(ri_arr)
                    sync_delay = max(0, sync_delay)

                    # Get throughput at join
                    # For open classes, use arrival rate from source or TN from parallel stations
                    if XN[r] > 0:
                        tput = XN[r]
                    else:
                        # For open networks, XN may be 0 - use TN from parallel stations
                        tput = 0.0
                        for station_idx in parallel_stations:
                            if TN[station_idx, r] > 0:
                                tput = TN[station_idx, r]
                                break
                        # If still 0, try to get arrival rate from source
                        if tput == 0:
                            source_stations = self._get_source_stations()
                            for src_idx in source_stations:
                                if self.rates[src_idx, r] > 0:
                                    tput = self.rates[src_idx, r]
                                    break

                    # Set join metrics
                    # RespT at join = synchronization delay
                    RN[join_station, r] = sync_delay
                    # QLen = tput * respT * num_branches (jobs from each branch waiting)
                    QN[join_station, r] = tput * sync_delay * num_branches
                    # Throughput at join = same as system throughput
                    TN[join_station, r] = tput
                    # Arrival rate at join = throughput * num_branches (jobs from all branches merge)
                    AN[join_station, r] = tput * num_branches
                    # Utilization stays 0 for join

                    # Update parallel station arrival rates
                    # Arrival rate = throughput / num_branches (visit ratio semantics)
                    for station_idx in parallel_stations:
                        AN[station_idx, r] = tput / num_branches

    def _find_parallel_branches(self, fork_idx: int, join_idx: int) -> List[int]:
        """Find stations on parallel branches between a fork and join.

        Uses the connection matrix to find nodes that are direct successors
        of the fork and eventually lead to the join.

        Returns a flat list of all station indices for backward compatibility.
        """
        branches = self._find_branch_structure(fork_idx, join_idx)
        # Flatten all branches into a single list
        parallel_stations = []
        for branch in branches:
            for station in branch:
                if station not in parallel_stations:
                    parallel_stations.append(station)
        return parallel_stations

    def _find_branch_structure(self, fork_idx: int, join_idx: int) -> List[List[int]]:
        """Find the branch structure between a fork and join.

        Returns a list of branches, where each branch is a list of station indices.
        For asymmetric fork-join, each branch may contain multiple serial stations.
        """
        from ..api.sn.network_struct import NodeType
        import numpy as np

        branches = []
        sn = self._sn

        # Use connection matrix to find direct successors of the fork
        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
            connmatrix = np.array(sn.connmatrix)

            # Get direct successors of the fork node
            fork_successors = np.where(connmatrix[fork_idx, :] > 0)[0]

            # For each successor, trace the path to see if it reaches the join
            for succ_node in fork_successors:
                # BFS to find all nodes on this branch until we hit the join
                visited = set()
                queue = [succ_node]
                branch_nodes = []

                while queue:
                    current = queue.pop(0)
                    if current in visited:
                        continue
                    visited.add(current)

                    if current == join_idx:
                        # Reached the join, this is a valid branch
                        break

                    branch_nodes.append(current)

                    # Add successors to queue
                    successors = np.where(connmatrix[current, :] > 0)[0]
                    for s in successors:
                        if s not in visited:
                            queue.append(s)

                # Convert branch nodes to station indices
                branch_stations = []
                for node in branch_nodes:
                    if node < len(sn.nodeToStation):
                        station = sn.nodeToStation[node]
                        if station >= 0 and station not in branch_stations:
                            # Check it's a queue/delay station (has service)
                            if station < len(self.station_types):
                                st = self.station_types[station]
                                if st is not None:
                                    st_val = st.value if hasattr(st, 'value') else int(st)
                                    if st_val in (NodeType.QUEUE, NodeType.DELAY):
                                        branch_stations.append(station)

                if branch_stations:
                    branches.append(branch_stations)
        else:
            # Fallback: treat all queues as single-station branches
            for i in range(self.nstations):
                if self.station_types and i < len(self.station_types):
                    st = self.station_types[i]
                    if st is not None:
                        st_val = st.value if hasattr(st, 'value') else int(st)
                        if st_val == NodeType.QUEUE:
                            branches.append([i])

        return branches

    def _compute_sync_delay(self, path_times: np.ndarray) -> float:
        """
        Compute synchronization delay using Heidelberger-Trivedi formula.

        For K parallel branches with response times r_1, ..., r_K,
        the expected maximum E[max(r_1, ..., r_K)] is computed using
        inclusion-exclusion with exponential approximation.

        Args:
            path_times: Array of response times for parallel paths

        Returns:
            Expected maximum (synchronization point) time
        """
        from itertools import combinations

        if len(path_times) == 0:
            return 0.0
        if len(path_times) == 1:
            return path_times[0]

        # Convert to rates (1/response_time)
        path_times = np.asarray(path_times)
        # Avoid division by zero
        path_times = np.maximum(path_times, 1e-10)
        lambdai = 1.0 / path_times

        d0 = 0.0
        parallel_branches = len(lambdai)

        for pow_val in range(parallel_branches):
            # Get all combinations of (pow_val + 1) elements
            for combo in combinations(range(parallel_branches), pow_val + 1):
                combo_sum = np.sum(lambdai[list(combo)])
                if combo_sum > 0:
                    d0 += ((-1) ** pow_val) * (1.0 / combo_sum)

        return d0

    def _find_node_type_indices(self, nodetype_list, target_type) -> np.ndarray:
        """
        Find indices of nodes with a specific type.

        Properly handles enum comparisons with lists of NodeType values.

        Args:
            nodetype_list: List or array of NodeType values
            target_type: The NodeType to search for

        Returns:
            NumPy array of indices where nodetype matches target_type
        """
        indices = []
        target_val = target_type.value if hasattr(target_type, 'value') else target_type
        for i, nt in enumerate(nodetype_list):
            nt_val = nt.value if hasattr(nt, 'value') else nt
            if nt_val == target_val:
                indices.append(i)
        return np.array(indices, dtype=int)

