"""Estimator implementations dispatched from ParamEstimator.

Each function takes (self, nodes) where self is a ParamEstimator instance.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""
import copy

import numpy as np
from scipy.optimize import nnls, minimize

from line_solver import (MetricType, SchedStrategy, Network, Delay, Queue, Exp,
                         ClosedClass, SolverMVA, SolverCTMC, EventType, Source)

from line_solver.inference.api.infer_qmle import infer_qmle
from line_solver.inference.api.infer_mlps import infer_mlps
from line_solver.inference.api.infer_fmlps import infer_fmlps
from line_solver.inference.api.infer_gibbs import infer_gibbs as _infer_gibbs
from line_solver.inference.api.infer_compute_ql_at_arrival import infer_compute_ql_at_arrival


def estimator_ubr(se, nodes):
    """Utilization-Based Regression."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()

    avg_aggr_util = None
    if np.isfinite(node.getNumberOfServers()):
        U = se.get_aggr_util(node)
        if U is not None:
            avg_aggr_util = U.data * node.getNumberOfServers()

    is_util_known = np.zeros(sn.nclasses, dtype=bool)
    avg_util = [None] * sn.nclasses
    avg_arvr = [None] * sn.nclasses
    for r in range(sn.nclasses):
        jc = se.model.getClasses()[r]
        Ur = se.get_util(node, jc)
        if Ur is not None:
            avg_util[r] = Ur.data * node.getNumberOfServers()
            is_util_known[r] = True
        arvr = se.get_arvr(node, jc)
        if arvr is None:
            raise ValueError(f'Arrival rate data for class {r + 1} is missing.')
        avg_arvr[r] = arvr.data

    # Check dimensions match
    n = len(avg_arvr[0])
    avg_a = np.column_stack(avg_arvr)
    sum_ur = np.zeros(n)
    for r in range(sn.nclasses):
        if is_util_known[r]:
            sum_ur += avg_util[r]

    est_val = np.zeros(sn.nclasses)

    # Known per-class utilization
    for r in np.where(is_util_known)[0]:
        x, _ = nnls(avg_a[:, r:r + 1], avg_util[r])
        est_val[r] = x[0]

    # Unknown classes: use aggregate utilization
    unknown = ~is_util_known
    if np.any(unknown) and avg_aggr_util is not None:
        residual_util = avg_aggr_util - sum_ur
        est_unknown, _ = nnls(avg_a[:, unknown], residual_util)
        est_val[unknown] = est_unknown

    return est_val


def estimator_ubo(se, nodes):
    """Utilization-Based Optimization (Liu et al. 2006)."""
    sn = se.model.refreshStruct() or se.model.getStruct()

    avg_u_list = []
    avg_arvr_dict = {}
    avg_respt_dict = {}

    for n_idx in range(len(nodes)):
        node = nodes[n_idx]
        if np.isfinite(node.getNumberOfServers()):
            U = se.get_aggr_util(node)
            if U is not None:
                avg_u_list.append(U.data * node.getNumberOfServers())

        for r in range(sn.nclasses):
            jc = se.model.getClasses()[r]
            arvr = se.get_arvr(node, jc)
            if arvr is None:
                raise ValueError(f'Arrival rate data for node {n_idx + 1} class {r + 1} missing.')
            avg_arvr_dict[(n_idx, r)] = arvr.data
            respt = se.get_respt(node, jc)
            if respt is None:
                raise ValueError(f'Response time data for node {n_idx + 1} class {r + 1} missing.')
            avg_respt_dict[(n_idx, r)] = respt.data

    M_nodes = len(nodes)
    R = sn.nclasses
    n_samples = len(avg_arvr_dict[(0, 0)])

    cpu_util = np.column_stack(avg_u_list) if avg_u_list else np.zeros((n_samples, M_nodes))
    avg_a = np.zeros((n_samples, M_nodes, R))
    avg_r = np.zeros((n_samples, M_nodes, R))
    for n_idx in range(M_nodes):
        for r in range(R):
            avg_a[:, n_idx, r] = avg_arvr_dict[(n_idx, r)]
            avg_r[:, n_idx, r] = avg_respt_dict[(n_idx, r)]

    return _ubo_qp(cpu_util, avg_r, avg_a)


def _ubo_qp(cpu_util, r_avg_times, avg_arvr):
    """Solve UBO as QP."""
    # Clean NaN
    valid = ~np.any(np.isnan(cpu_util), axis=1)
    cpu_util = cpu_util[valid]
    r_avg_times = r_avg_times[valid]
    avg_arvr = avg_arvr[valid]

    # Remove zero throughput
    valid2 = np.sum(np.sum(avg_arvr, axis=2), axis=1) > 0
    cpu_util = cpu_util[valid2]
    r_avg_times = r_avg_times[valid2]
    avg_arvr = avg_arvr[valid2]

    N_exp = cpu_util.shape[0]
    M = cpu_util.shape[1]
    R = r_avg_times.shape[2]
    MR = M * R

    H = np.zeros((MR, MR))
    h = np.zeros(MR)

    for n in range(N_exp):
        rho_n = cpu_util[n, :]
        beta_n = 1.0 / (1.0 - rho_n + 1e-10)

        lambda_n = avg_arvr[n, :, :]  # M x R
        R_n = r_avg_times[n, :, :]  # M x R

        lambda_r = np.sum(lambda_n, axis=0)  # R
        w_n = lambda_r / (np.sum(lambda_r) + 1e-10)

        E_n = np.sum(R_n, axis=0)  # R

        # A_delta: R x MR
        A_delta = np.zeros((R, MR))
        for r in range(R):
            A_delta[r, r * M:(r + 1) * M] = beta_n

        # A_eps: M x MR
        A_eps = np.zeros((M, MR))
        for i in range(M):
            for r in range(R):
                A_eps[i, r * M + i] = lambda_n[i, r]

        W_n = np.diag(w_n)

        H += 2 * (A_delta.T @ W_n @ A_delta + A_eps.T @ A_eps)
        h -= 2 * (A_delta.T @ W_n @ E_n + A_eps.T @ rho_n)

    # Solve QP: min 0.5 * s^T H s + h^T s, s.t. s >= 0
    # Use scipy minimize with bounds
    x0 = np.ones(MR) * 0.1
    bounds = [(0, None)] * MR

    def qp_obj(s):
        return 0.5 * s @ H @ s + h @ s

    def qp_grad(s):
        return H @ s + h

    result = minimize(qp_obj, x0, jac=qp_grad, method='L-BFGS-B', bounds=bounds)
    return result.x.reshape(M, R, order='F')


def estimator_erps(se, nodes):
    """Extended Regression for Processor Sharing."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()

    if node.getSchedStrategy() != SchedStrategy.PS:
        raise ValueError('ERPS is only available for processor sharing stations.')

    R = sn.nclasses
    avg_respt = [None] * R
    avg_aqlen = [None] * R

    for r in range(R):
        jc = se.model.getClasses()[r]

        respt = se.get_respt(node, jc)
        if respt is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        avg_respt[r] = respt.data

        # Get aggregate queue-length conditional on class-r arrivals
        try:
            from line_solver import Event, EventType as ET
            arv_event = Event(ET.ARV, node, jc)
            aqlen = se.get_aggr_qlen(node, arv_event)
        except (ImportError, TypeError):
            aqlen = se.get_aggr_qlen(node)

        if aqlen is None:
            raise ValueError(f'Arrival queue-length data for class {r + 1} missing.')
        avg_aqlen[r] = aqlen.data

    # Estimate average busy cores
    busy_cores = np.concatenate([np.sum(avg_aqlen[r].reshape(-1, 1), axis=1) for r in range(R)])
    avg_busy_cores = min(np.mean(busy_cores), node.getNumberOfServers())

    est_val = np.zeros(R)
    for r in range(R):
        resp_times = avg_respt[r]
        total_ql = np.sum(avg_aqlen[r].reshape(-1, 1), axis=1) / avg_busy_cores
        x, _ = nnls(total_ql.reshape(-1, 1), resp_times)
        est_val[r] = x[0]

    return est_val


def estimator_ekf(se, nodes):
    """Extended Kalman Filter resource demand estimator."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()

    avg_u = None
    if np.isfinite(node.getNumberOfServers()):
        U = se.get_aggr_util(node)
        if U is not None:
            avg_u = U.data * node.getNumberOfServers()

    avg_arvr = []
    avg_respt = []
    for r in range(sn.nclasses):
        jc = se.model.getClasses()[r]
        arvr = se.get_arvr(node, jc)
        if arvr is None:
            raise ValueError(f'Arrival rate data for class {r + 1} missing.')
        avg_arvr.append(arvr.data)
        respt = se.get_respt(node, jc)
        if respt is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        avg_respt.append(respt.data)

    avg_a = np.column_stack(avg_arvr)
    avg_r = np.column_stack(avg_respt)

    solver_type = se.options.get('solver', SolverMVA)
    x0 = se.options.get('x0', None)

    return _ekf_data(avg_u, avg_r, avg_a, node.getNumberOfServers(),
                     se.options['iter_max'], se.model, node, solver_type, x0)


def _ekf_data(cpu_util, r_avg_times, avg_arvr, num_servers, iter_max,
              model, node, solver_type, x0_init):
    # Clean NaN
    valid = ~np.isnan(cpu_util)
    cpu_util = cpu_util[valid]
    r_avg_times = r_avg_times[valid]
    avg_arvr = avg_arvr[valid]

    # Remove zero throughput
    valid2 = np.sum(avg_arvr, axis=1) > 0
    cpu_util = cpu_util[valid2]
    r_avg_times = r_avg_times[valid2]
    avg_arvr = avg_arvr[valid2]

    R = r_avg_times.shape[1]

    if x0_init is not None:
        x = np.asarray(x0_init, dtype=float).flatten()
    else:
        x = np.random.rand(R) * np.max(r_avg_times, axis=0)

    p = np.diag(x ** 2)
    m_cov_noise = np.diag(np.full(R + 1, 0.01))
    p_cov_noise = np.eye(R) * 0.001

    # Initialize solver
    sn = model.refreshStruct() or model.getStruct()
    st_idx = node.getStationIndex() - 1  # 0-based

    solver = solver_type(model)
    solver_analyzer = _get_solver_analyzer(solver, model, sn, st_idx)

    N_exp = cpu_util.shape[0]
    step_bound = 0.6
    a_min = np.zeros(R)
    a_max = np.full(R, np.inf)

    for n in range(min(N_exp, iter_max)):
        x_n = x.copy()
        P_n = p + p_cov_noise

        step_util = cpu_util[n]
        step_response = r_avg_times[n, :]

        z_n, sn = _get_predicted_measurement(x_n, R, sn, st_idx, solver_analyzer)
        H_n = _get_jacobian(x_n, R, sn, st_idx, solver_analyzer, z_n)
        z = np.concatenate([step_response, [step_util]])
        y_n = z - z_n

        H_nT = H_n.T
        S_n = H_n @ P_n @ H_nT + m_cov_noise
        K_n = P_n @ H_nT @ np.linalg.inv(S_n)

        x = x_n + K_n @ y_n
        x_lower = step_bound * a_min + (1 - step_bound) * x
        x_upper = step_bound * a_max + (1 - step_bound) * x
        x = np.minimum(x_upper, np.maximum(x_lower, x))
        if np.sum(x) < 0:
            x = -x
        p = (np.eye(R) - K_n @ H_n) @ P_n

    return x


def _get_solver_analyzer(solver, model, sn, st_idx):
    """Return a callable that takes sn and returns (Q, U, R, T)."""
    from line_solver.api.solvers.mva.analyzers import solver_mva_analyzer

    def analyzer(sn_arg):
        try:
            result = solver_mva_analyzer(sn_arg)
            return result.QN, result.UN, result.RN, result.TN
        except Exception:
            solver_obj = SolverMVA(model)
            solver_obj.runAnalyzer()
            Q = solver_obj.result.QN
            U = solver_obj.result.UN
            R = solver_obj.result.RN
            T = solver_obj.result.TN
            return Q, U, R, T

    return analyzer


def _get_predicted_measurement(x, R, sn, st_idx, solver_analyzer):
    from line_solver.api.sn.transforms import sn_set_service
    for c in range(R):
        sn_set_service(sn, st_idx, c, 1.0 / x[c], 1.0)
    Q, U, RN, T = solver_analyzer(sn)
    h = np.zeros(R + 1)
    for c in range(R):
        h[c] = RN[st_idx, c]
    h[R] = np.sum(U[st_idx, :])
    return h, sn


def _get_jacobian(x, R, sn, st_idx, solver_analyzer, h0):
    from line_solver.api.sn.transforms import sn_set_service
    delta = 1e-6
    Hx = np.zeros((R + 1, R))
    for c in range(R):
        x_pert = x.copy()
        x_pert[c] += delta
        sn_pert = copy.deepcopy(sn)
        for cc in range(R):
            sn_set_service(sn_pert, st_idx, cc, 1.0 / x_pert[cc], 1.0)
        Q, U_pert, R_pert, T = solver_analyzer(sn_pert)
        h_pert = np.zeros(R + 1)
        for cc in range(R):
            h_pert[cc] = R_pert[st_idx, cc]
        h_pert[R] = np.sum(U_pert[st_idx, :])
        Hx[:, c] = (h_pert - h0) / delta
    return Hx


def estimator_mcmc(se, nodes):
    """Gibbs Sampling MCMC-based optimization."""
    sn = se.model.refreshStruct() or se.model.getStruct()

    Nopen = se.options.get('openPopulation', 100)

    P_pop = np.zeros(sn.nclasses)
    for r in range(sn.nclasses):
        if sn.njobs[r] < np.inf:
            P_pop[r] = sn.njobs[r]
        else:
            P_pop[r] = Nopen

    Z = np.zeros(sn.nclasses)
    classes = se.model.getClasses()
    all_nodes = se.model.getNodes()
    for nd in all_nodes:
        if isinstance(nd, Delay):
            for r in range(sn.nclasses):
                if sn.njobs[r] < np.inf:
                    Z[r] += nd.getService(classes[r]).getMean()
        elif isinstance(nd, Source):
            for r in range(sn.nclasses):
                if sn.njobs[r] == np.inf:
                    lambda_r = 1.0 / nd.getService(classes[r]).getMean()
                    Z[r] = P_pop[r] / lambda_r

    # Collect aggregate queue-length data
    avg_ql_list = []
    for n_idx in range(len(nodes)):
        node = nodes[n_idx]
        aqlen = se.get_aggr_qlen(node)
        if aqlen is None:
            raise ValueError(f'Transient queue-length data for node {n_idx + 1} missing.')
        avg_ql_list.append(aqlen.data)

    avg_ql = np.column_stack(avg_ql_list)
    experiments = avg_ql.shape[0]
    avg_ql_mean = np.mean(avg_ql, axis=0).reshape(len(nodes), -1)

    return _mcmc_data(avg_ql_mean, sn.visits, experiments, se.options['iter_max'], P_pop, Z)


def _mcmc_data(avg_ql, visits, experiments, iter_max, P_pop, Z):
    from line_solver.api.pfqn.mva import pfqn_bs
    from line_solver.api.pfqn.asymptotic import pfqn_mci

    M = avg_ql.shape[0]
    R = avg_ql.shape[1] if avg_ql.ndim > 1 else 1
    if avg_ql.ndim == 1:
        avg_ql = avg_ql.reshape(M, 1)

    S = 100
    integral_range = [0, np.max(avg_ql)]
    theta_step = integral_range[1] / 400

    theta = np.zeros((S + 1, M, R))
    steps = np.arange(integral_range[0], integral_range[1] + theta_step, theta_step)

    for s in range(S):
        sample_theta = theta[s].copy()
        for i in range(M):
            for c in range(R):
                log_posteriors = np.zeros(len(steps))
                log_prior = np.log(theta_step / (integral_range[1] - integral_range[0] + 1e-300))
                try:
                    g_nc, _, _ = pfqn_mci(sample_theta, P_pop, Z, experiments, variant='imci')
                except Exception:
                    g_nc = 1.0

                for st in range(len(steps)):
                    log_posteriors[st] = (experiments * avg_ql[i, c] * np.log(steps[st] + 1e-300)
                                          - experiments * np.log(g_nc + 1e-300) + log_prior)

                probs = np.exp(log_posteriors - np.max(log_posteriors))
                probs /= np.sum(probs)

                cum_prob = np.cumsum(probs)
                u = np.random.rand()
                index = np.searchsorted(cum_prob, u)
                if index >= len(steps):
                    index = len(steps) - 1
                sample_theta[i, c] = steps[index]

        theta[s + 1] = sample_theta

    # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
    visit_per_class = np.ones((M, R))
    for r in range(R):
        if r in visits:
            v = visits[r]
            if v.ndim == 2 and v.shape[0] >= M:
                visit_per_class[:, r] = v[:M, r]

    # Avoid division by zero
    visit_per_class = np.where(visit_per_class > 0, visit_per_class, 1.0)

    # Divide all theta samples by visit ratios
    for s in range(S + 1):
        theta[s] = theta[s] / visit_per_class

    cutoff = S // 2
    theta_avg = np.mean(theta[cutoff:], axis=0)

    return theta_avg


def estimator_mle(se, nodes):
    """Maximum Likelihood Estimation."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()

    avg_u = None
    if np.isfinite(node.getNumberOfServers()):
        U = se.get_aggr_util(node)
        if U is not None:
            avg_u = U.data * node.getNumberOfServers()

    avg_arvr = []
    avg_respt = []
    for r in range(sn.nclasses):
        jc = se.model.getClasses()[r]
        arvr = se.get_arvr(node, jc)
        if arvr is None:
            raise ValueError(f'Arrival rate data for class {r + 1} missing.')
        avg_arvr.append(arvr.data)
        respt = se.get_respt(node, jc)
        if respt is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        avg_respt.append(respt.data)

    avg_a = np.column_stack(avg_arvr)
    avg_r = np.column_stack(avg_respt)

    solver_type = se.options.get('solver', SolverMVA)

    return _mle_data(avg_u, avg_r, avg_a, node.getNumberOfServers(),
                     se.options['iter_max'], se.model, node, solver_type)


def _mle_data(cpu_util, r_avg_times, avg_arvr, num_servers, iter_max,
              model, node, solver_type):
    valid = ~np.isnan(cpu_util)
    cpu_util = cpu_util[valid]
    r_avg_times = r_avg_times[valid]
    avg_arvr = avg_arvr[valid]

    valid2 = np.sum(avg_arvr, axis=1) > 0
    cpu_util = cpu_util[valid2]
    r_avg_times = r_avg_times[valid2]
    avg_arvr = avg_arvr[valid2]

    R = r_avg_times.shape[1]

    sn = model.refreshStruct() or model.getStruct()
    st_idx = node.getStationIndex() - 1

    solver = solver_type(model)
    solver_analyzer = _get_solver_analyzer(solver, model, sn, st_idx)

    x0 = np.random.rand(R) * np.max(r_avg_times, axis=0)

    N_exp = cpu_util.shape[0]
    w = avg_arvr / (np.sum(avg_arvr, axis=1, keepdims=True) + 1e-10)

    def objfun(x):
        from line_solver.api.sn.transforms import sn_set_service
        for c in range(R):
            sn_set_service(sn, st_idx, c, 1.0 / x[c], 1.0)
        Q, U_pred, R_pred, T = solver_analyzer(sn)
        pred_r = R_pred[st_idx, :R]
        pred_u = np.sum(U_pred[st_idx, :R])

        delta_j = np.tile(pred_r, (N_exp, 1)) - r_avg_times
        epsi = pred_u * np.ones(N_exp) - cpu_util
        f = np.sum(w * delta_j ** 2) + np.sum(epsi ** 2)
        return f

    x_lb = np.full(R, 1e-8)
    x_ub = np.max(r_avg_times, axis=0)
    bounds = list(zip(x_lb, x_ub))

    result = minimize(objfun, x0, method='L-BFGS-B', bounds=bounds,
                      options={'maxiter': iter_max})
    return result.x


def estimator_rnn(se, nodes):
    """Explainable RNN estimation."""
    import torch
    import torch.nn as nn
    from line_solver.inference.lang.rnn_layer import QueueNetworkLearningRNNLayer

    sn = se.model.refreshStruct() or se.model.getStruct()
    all_nodes = se.model.getNodes()

    ql_ts = {}
    ql_trace = {}
    num_servers = []
    min_sample_count = float('inf')

    for n_idx, nd in enumerate(all_nodes):
        num_servers.append(nd.getNumberOfServers())
        for r in range(sn.nclasses):
            jc = se.model.getClasses()[r]
            samples = se.get_qlen(nd, jc)
            if samples is None:
                raise ValueError(f'Queue-length data for node {n_idx + 1} class {r + 1} missing.')
            if not isinstance(samples, list):
                samples = [samples]
            for d, ql_data in enumerate(samples):
                if d not in ql_ts:
                    ql_ts[d] = {}
                    ql_trace[d] = {}
                ql_ts[d][(n_idx, r)] = ql_data.t
                ql_trace[d][(n_idx, r)] = ql_data.data
                if len(ql_data.data) < min_sample_count:
                    min_sample_count = len(ql_data.data)

    M = len(all_nodes)
    R = sn.nclasses

    traces_list = []
    for d in sorted(ql_ts.keys()):
        trace = np.zeros((int(min_sample_count), M, R + 1))
        for n_idx in range(M):
            for r in range(R):
                data = ql_trace[d][(n_idx, r)][:int(min_sample_count)]
                ts = ql_ts[d][(n_idx, r)][:int(min_sample_count)]
                trace[:, n_idx, 0] = ts
                trace[:, n_idx, r + 1] = data
        traces_list.append(trace)

    traces = np.stack(traces_list, axis=0)

    return _rnn_data(traces, np.array(num_servers))


def _rnn_data(avg_ql, num_servers):
    import torch
    import torch.nn as nn
    from line_solver.inference.lang.rnn_layer import QueueNetworkLearningRNNLayer

    M = avg_ql.shape[2]
    R = avg_ql.shape[3] - 1
    S = avg_ql.shape[1]
    trace_count = avg_ql.shape[0]

    num_epochs = 2
    num_iters_per_epoch = 50

    layer = QueueNetworkLearningRNNLayer(M, R, num_servers)
    optimizer = torch.optim.Adam(layer.parameters(), lr=0.1)

    for epoch in range(num_epochs):
        for i in range(num_iters_per_epoch):
            layer.reset_state()
            exp_idx = np.random.randint(trace_count)
            trace = avg_ql[exp_idx]

            X = torch.tensor(trace, dtype=torch.float32)
            T = torch.tensor(trace, dtype=torch.float32)

            optimizer.zero_grad()
            Y = layer(X)

            # Loss: max absolute percentage error
            pred_err = torch.abs(T[:, :, 1:] - Y[:, :, 1:])
            N_mean = torch.mean(torch.sum(X[:, :, 1:], dim=1))
            max_err = torch.max(torch.sum(pred_err, dim=0) / (2.0 * N_mean + 1e-10))
            loss = 100 * max_err

            loss.backward()
            optimizer.step()

    mu = torch.abs(layer.mu).detach().numpy()
    return 1.0 / mu.flatten()


def estimator_mlps(se, nodes):
    """Maximum Likelihood for Processor Sharing."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()
    R = sn.nclasses

    if node.getSchedStrategy() != SchedStrategy.PS:
        raise ValueError('MLPS is only available for processor sharing stations.')

    has_open = np.any(sn.njobs == np.inf)
    if has_open:
        eq_model, eq_node = se._build_closed_equivalent_for_ps(node)
    else:
        eq_model = se.model
        eq_node = node

    rt_all, class_all, at_all = [], [], []
    for r in range(R):
        jc = se.model.getClasses()[r]
        arv_data = se.get_arvr(node, jc)
        if arv_data is None:
            raise ValueError(f'Arrival timestamp data for class {r + 1} missing.')
        if not arv_data.is_trace():
            raise ValueError('MLPS requires trace-format arrival data.')

        rt_data = se.get_respt(node, jc)
        if rt_data is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        if not rt_data.is_trace():
            raise ValueError('MLPS requires trace-format response time data.')

        n_samples = len(rt_data.data)
        rt_all.append(rt_data.data)
        at_all.append(arv_data.data)
        class_all.append(np.full(n_samples, r + 1, dtype=int))

    rt_arr = np.concatenate(rt_all)
    at_arr = np.concatenate(at_all)
    class_arr = np.concatenate(class_all)

    n = len(at_arr)
    jobid = np.arange(1, n + 1)
    ql = infer_compute_ql_at_arrival(at_arr, jobid, rt_arr, jobid, class_arr, R)

    sort_idx = np.argsort(at_arr)
    rt_sorted = rt_arr[sort_idx]
    class_sorted = class_arr[sort_idx]
    ql = ql[sort_idx, :]

    valid = rt_sorted > 0
    rt_sorted = rt_sorted[valid]
    class_sorted = class_sorted[valid]
    ql = ql[valid, :]

    return infer_mlps(eq_model, eq_node, rt_sorted, class_sorted, ql)


def estimator_fmlps(se, nodes):
    """Fluid Maximum Likelihood for Processor Sharing."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()
    R = sn.nclasses

    if node.getSchedStrategy() != SchedStrategy.PS:
        raise ValueError('FMLPS is only available for processor sharing stations.')

    has_open = np.any(sn.njobs == np.inf)
    if has_open:
        eq_model, eq_node = se._build_closed_equivalent_for_ps(node)
        eq_sn = eq_model.refreshStruct() or eq_model.getStruct()
        W = int(np.sum(eq_sn.njobs))
    else:
        eq_model = se.model
        eq_node = node
        W = int(np.sum(sn.njobs))

    rt_all, class_all, at_all = [], [], []
    for r in range(R):
        jc = se.model.getClasses()[r]
        arv_data = se.get_arvr(node, jc)
        if arv_data is None:
            raise ValueError(f'Arrival timestamp data for class {r + 1} missing.')
        if not arv_data.is_trace():
            raise ValueError('FMLPS requires trace-format arrival data.')

        rt_data = se.get_respt(node, jc)
        if rt_data is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        if not rt_data.is_trace():
            raise ValueError('FMLPS requires trace-format response time data.')

        n_samples = len(rt_data.data)
        rt_all.append(rt_data.data)
        at_all.append(arv_data.data)
        class_all.append(np.full(n_samples, r + 1, dtype=int))

    rt_arr = np.concatenate(rt_all)
    at_arr = np.concatenate(at_all)
    class_arr = np.concatenate(class_all)

    n = len(at_arr)
    jobid = np.arange(1, n + 1)
    ql = infer_compute_ql_at_arrival(at_arr, jobid, rt_arr, jobid, class_arr, R)

    sort_idx = np.argsort(at_arr)
    rt_sorted = rt_arr[sort_idx]
    class_sorted = class_arr[sort_idx]
    ql = ql[sort_idx, :]

    valid = rt_sorted > 0
    rt_sorted = rt_sorted[valid]
    class_sorted = class_sorted[valid]
    ql = ql[valid, :]

    return infer_fmlps(eq_model, eq_node, rt_sorted, class_sorted, ql, W)


def estimator_qmle(se, nodes):
    """Quick Maximum Likelihood Estimation."""
    sn = se.model.refreshStruct() or se.model.getStruct()
    R = sn.nclasses
    M = len(nodes)

    Nopen = se.options.get('openPopulation', 100)

    N = np.zeros(R)
    for r in range(R):
        if sn.njobs[r] < np.inf:
            N[r] = sn.njobs[r]
        else:
            N[r] = Nopen

    Z = np.zeros(R)
    classes = se.model.getClasses()
    all_nodes = se.model.getNodes()
    for nd in all_nodes:
        if isinstance(nd, Delay):
            for r in range(R):
                if sn.njobs[r] < np.inf:
                    Z[r] += nd.getService(classes[r]).getMean()
        elif isinstance(nd, Source):
            for r in range(R):
                if sn.njobs[r] == np.inf:
                    lambda_r = 1.0 / nd.getService(classes[r]).getMean()
                    Z[r] = N[r] / lambda_r

    Q = np.zeros((M, R))
    for n_idx in range(M):
        for r in range(R):
            jc = classes[r]
            ql_data = se.get_qlen(nodes[n_idx], jc)
            if ql_data is None:
                raise ValueError(f'Queue-length data for node {n_idx + 1} class {r + 1} missing.')
            if isinstance(ql_data, list):
                ql_data = ql_data[0]
            Q[n_idx, r] = np.mean(ql_data.data)

    return infer_qmle(Q, N, Z)


def estimator_gibbs(se, nodes):
    """Gibbs Sampling from trace data."""
    node = nodes[0]
    sn = se.model.refreshStruct() or se.model.getStruct()
    R = sn.nclasses
    nb_cores = node.getNumberOfServers()

    # Build legacy data format
    data = [[] for _ in range(6)]
    for _ in range(R):
        for j in range(6):
            data[j].append(None)

    for r in range(R):
        jc = se.model.getClasses()[r]

        arv_data = se.get_arvr(node, jc)
        if arv_data is None:
            raise ValueError(f'Arrival data for class {r + 1} missing.')
        if arv_data.is_trace():
            data[2][r] = arv_data.data * 1000  # s -> ms
        else:
            raise ValueError('Gibbs requires trace-format arrival data.')

        rt_data = se.get_respt(node, jc)
        if rt_data is None:
            raise ValueError(f'Response time data for class {r + 1} missing.')
        if rt_data.is_trace():
            data[3][r] = rt_data.data
        else:
            raise ValueError('Gibbs requires trace-format response time data.')

        tput_data = se.get_tput(node, jc)
        if tput_data is None:
            raise ValueError(f'Throughput data for class {r + 1} missing.')
        data[5][r] = tput_data.data

    return _infer_gibbs(data, nb_cores, se.options['tol'])
