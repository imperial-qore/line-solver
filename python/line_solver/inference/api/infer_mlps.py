import numpy as np
from scipy.optimize import minimize
from scipy.linalg import expm

from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy
from line_solver import SolverCTMC, EventType
from line_solver.api.sn.transforms import sn_set_service
from line_solver.api.solvers.ctmc import solver_ctmc, SolverCTMCOptions


def infer_mlps(model, node, rt, cls, ql):
    """MLPS demand estimation using sn struct-level operations.

    Estimates service demands at a PS queue using Maximum Likelihood
    for Processor Sharing. Pre-builds augmented CTMC models for each
    unique (tagClass, aQueue) combination, then uses sn_set_service
    + solver_ctmc directly inside the optimization loop.

    Args:
        model: LINE Network model
        node: PS queue node
        rt: response time samples (n,)
        cls: class of each sample (n,), 1-based
        ql: queue lengths at arrival (n x R)

    Returns:
        demand_est: 1-D array of estimated demands (R,)

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    rt = np.asarray(rt, dtype=float).flatten()
    cls = np.asarray(cls, dtype=int).flatten()
    ql = np.asarray(ql, dtype=float)

    sn = model.refreshStruct() or model.getStruct()
    R = sn.nclasses
    n_cores = node.getNumberOfServers()

    # Get delay rates from model
    delay_idx = -1
    for i in range(sn.nstations):
        if sn.sched[i] == SchedStrategy.INF or sn.sched[i] == SchedStrategy.INF.value:
            delay_idx = i
            break
    mu_z = np.zeros(R)
    for k in range(R):
        mu_z[k] = sn.rates[delay_idx, k]

    # Initial point estimate
    mean_ql = np.mean(np.sum(ql, axis=1))
    v_tilde = min(mean_ql, n_cores)
    x0 = np.zeros(R)
    for j in range(R):
        mask = cls == (j + 1)
        if np.any(mask):
            x0[j] = v_tilde * np.mean(rt[mask]) / mean_ql
        else:
            x0[j] = 1e-3

    x_lb = np.full(R, 1e-8)
    x_ub = np.full(R, np.max(rt))

    new_r = R + 1
    aug_mu_z = np.append(mu_z, 0.0)

    # Pre-build augmented CTMC models for all unique (tagClass, aQueue)
    unique_tc = np.unique(cls)
    unique_ql = np.unique(ql, axis=0)

    ctmc_opts = SolverCTMCOptions(gen_method='sync')

    prebuilt = {}
    for tc in unique_tc:
        aug_mu_z[new_r - 1] = mu_z[tc - 1]
        for v_idx in range(unique_ql.shape[0]):
            aq = unique_ql[v_idx, :]

            # Population: move 1 job from tagClass to tagged class
            N_pop = aq.copy()
            N_pop[tc - 1] -= 1
            N_pop = np.append(N_pop, 1)

            # Build augmented Network model (once per combo)
            aug_model = Network('mlps_aug')
            aug_delay = Delay(aug_model, 'Think')
            aug_queue = Queue(aug_model, 'Queue1', SchedStrategy.PS)
            aug_queue.setNumberOfServers(n_cores)
            aug_classes = []
            for r in range(new_r):
                ac = ClosedClass(aug_model, f'Class{r + 1}', int(N_pop[r]), aug_delay, 0)
                aug_classes.append(ac)
                aug_delay.setService(ac, Exp(aug_mu_z[r]))
                aug_queue.setService(ac, Exp(1.0))  # placeholder
            P = aug_model.initRoutingMatrix()
            for r in range(new_r):
                P.set(aug_classes[r], aug_classes[r], aug_delay, aug_queue, 1.0)
                P.set(aug_classes[r], aug_classes[r], aug_queue, aug_delay, 1.0)
            aug_model.link(P)

            # Use SolverCTMC for state space, generator, and event filtering
            solver = SolverCTMC(aug_model)
            solver.runAnalyzer()
            inf_gen, event_filt = solver.getGenerator()
            ss_aggr = solver.getStateSpaceAggr()

            # Find tagged departure event indices
            queue_node_idx = aug_model.get_node_index(aug_queue)
            tagged_dep_idx = []
            if hasattr(solver, '_result') and hasattr(solver._result, 'event'):
                events = solver._result.event
                for e_idx, ev in enumerate(events):
                    if (hasattr(ev, 'active') and ev.active
                            and ev.active[0].node == queue_node_idx
                            and ev.active[0].jobclass == new_r
                            and ev.active[0].event == EventType.DEP):
                        tagged_dep_idx.append(e_idx)

            # Find subset where tagged job is at Queue
            queue_st_idx = aug_queue.getStationIndex() - 1  # 0-based
            tagged_col = queue_st_idx * new_r + (new_r - 1)
            subset = np.where(ss_aggr[:, tagged_col] == 1)[0]

            # Extract queue state space for subset
            queue_cols = slice(queue_st_idx * new_r, queue_st_idx * new_r + new_r)
            ss_queue = ss_aggr[subset][:, queue_cols]

            # Store the augmented sn struct for rate updates in objfun
            aug_sn = aug_model.refreshStruct() or aug_model.getStruct()
            queue_st_idx_0 = aug_queue.getStationIndex() - 1  # 0-based for sn_set_service

            cache_key = (int(tc), tuple(aq.astype(int)))
            prebuilt[cache_key] = {
                'sn': aug_sn,
                'queue_st_idx': queue_st_idx_0,
                'tagged_dep_idx': tagged_dep_idx,
                'subset': subset,
                'ss_queue': ss_queue,
                'N': N_pop,
                'tag_class': int(tc),
            }

    def objfun(x):
        TOL = 1e-6
        rates = 1.0 / x

        cache = {}
        for key, pb in prebuilt.items():
            tc_val = pb['tag_class']
            aug_rates = np.append(rates, rates[tc_val - 1])

            # Update service rates in cached sn struct
            sn_upd = pb['sn']
            for cc in range(new_r):
                sn_set_service(sn_upd, pb['queue_st_idx'], cc, aug_rates[cc])

            # Re-solve CTMC with updated rates
            ctmc_result = solver_ctmc(sn_upd, ctmc_opts)
            inf_gen = ctmc_result.infgen
            event_filt = ctmc_result.eventFilt

            # Build D1 from tagged departure event filters
            n_states = inf_gen.shape[0]
            D1 = np.zeros((n_states, n_states))
            if event_filt is not None:
                for di in pb['tagged_dep_idx']:
                    if di < len(event_filt):
                        D1 = D1 + event_filt[di]

            # Extract sub-generator: MAPQ1 = infGen - D1
            MAPQ1 = inf_gen - D1
            subset = pb['subset']
            A_sub = MAPQ1[np.ix_(subset, subset)]

            cache[key] = {
                'A': A_sub,
                'ss_queue': pb['ss_queue'],
                'N': pb['N'],
            }

        # Compute likelihoods
        ftemp = np.zeros(len(rt))
        for i in range(len(rt)):
            cache_key = (int(cls[i]), tuple(ql[i, :].astype(int)))
            cached = cache[cache_key]
            ftemp[i] = np.log(TOL + _eval_mlps_likelihood(
                cached['A'], cached['ss_queue'], cached['N'], rt[i]))
        return -np.sum(ftemp)

    bounds = list(zip(x_lb, x_ub))
    result = minimize(objfun, x0, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': 10000, 'ftol': 1e-10})
    return result.x


def _eval_mlps_likelihood(A, ss_queue, N, r_sam):
    """Compute MLPS likelihood from pre-built CTMC components."""
    n_states = A.shape[0]
    pie = np.zeros(n_states)

    # Find matching row
    N_arr = np.asarray(N).flatten()
    for i in range(ss_queue.shape[0]):
        if np.array_equal(ss_queue[i, :len(N_arr)], N_arr):
            pie[i] = 1.0
            break

    if np.sum(pie) == 0:
        return 0.0

    # MAP PDF: f(t) = pie * exp(A*t) * (-A*1)
    # where 1 is the column vector of ones
    try:
        exp_at = expm(A * r_sam)
        exit_rate = -A @ np.ones(n_states)
        like = pie @ exp_at @ exit_rate
        return max(0.0, float(like))
    except Exception:
        return 0.0
