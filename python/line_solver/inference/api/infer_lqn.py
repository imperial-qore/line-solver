# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.
"""LQN parameter identification via an Extended Kalman Filter.

Estimates hidden Layered Queueing Network parameters (activity host demands,
task think times) from measured performance data (response times, utilizations,
throughputs), following Zheng, Yang, Woodside, Litoiu, Iszlai, "Tracking
Time-Varying Parameters in Software Systems with Extended Kalman Filters",
CASCON 2005 (equations 1-9).

The parameter is modelled as a zero-mean random walk a_k = a_{k-1} + w and the
measurement as z_k = h(a_k) + v, where h is the LQN performance model. This is
the native-Python counterpart of the MATLAB ``infer_lqn`` family and the JAR
``jline.api.infer.InferLqn``; results match across codebases.

A parameter is specified as a dict {'type': 'hostdem'|'think', 'name': <name>};
an observation as {'metric': 'RespT'|'Util'|'Tput'|'QLen', 'name': <name>}.
"""
import numpy as np

__all__ = [
    'infer_lqn',
    'infer_lqn_ekf',
    'infer_lqn_jacobian',
    'infer_lqn_setparams',
    'infer_lqn_getparams',
    'infer_lqn_getobs',
    'infer_lqn_findbyname',
]


def infer_lqn_findbyname(container, name):
    """Return the first element of container whose .name equals name, else None."""
    for el in container:
        if getattr(el, 'name', None) == name:
            return el
    return None


def infer_lqn_setparams(model, param_spec, a):
    """Apply a parameter vector to a LayeredNetwork in place.

    param_spec[i] is a dict with keys 'type' ('hostdem' | 'think') and 'name'
    (activity name for 'hostdem', task name for 'think'). Means are injected as
    exponential distributions (SCV = 1).
    """
    a = np.asarray(a, dtype=float).ravel()
    if a.size != len(param_spec):
        raise ValueError('Length of parameter vector does not match param_spec.')
    for i, spec in enumerate(param_spec):
        val = float(a[i])
        typ = str(spec['type']).lower()
        if typ == 'hostdem':
            act = infer_lqn_findbyname(model.activities, spec['name'])
            if act is None:
                raise ValueError("Activity '%s' not found." % spec['name'])
            act.setHostDemand(val)
        elif typ == 'think':
            tsk = infer_lqn_findbyname(model.tasks, spec['name'])
            if tsk is None:
                raise ValueError("Task '%s' not found." % spec['name'])
            tsk.setThinkTime(val)
        else:
            raise ValueError("Unknown parameter type '%s'." % spec['type'])
    return model


def infer_lqn_getparams(model, param_spec):
    """Read the current values of the parameters named in param_spec."""
    a0 = np.zeros(len(param_spec))
    for i, spec in enumerate(param_spec):
        typ = str(spec['type']).lower()
        if typ == 'hostdem':
            act = infer_lqn_findbyname(model.activities, spec['name'])
            if act is None:
                raise ValueError("Activity '%s' not found." % spec['name'])
            a0[i] = act.getHostDemandMean()
        elif typ == 'think':
            tsk = infer_lqn_findbyname(model.tasks, spec['name'])
            if tsk is None:
                raise ValueError("Task '%s' not found." % spec['name'])
            a0[i] = tsk.getThinkTimeMean()
        else:
            raise ValueError("Unknown parameter type '%s'." % spec['type'])
    return a0


def infer_lqn_getobs(names, metrics, obs_spec):
    """Extract the observation vector selected by obs_spec.

    names is the element name list aligned with the metric vectors (as in
    LayeredNetworkStruct.names). metrics is a dict with keys QLen/Util/RespT/Tput
    (each a vector indexed by element). obs_spec[i] is a dict with keys 'metric'
    ('RespT'|'Util'|'Tput'|'QLen') and 'name' (element name).
    """
    names_list = list(names)
    z = np.zeros(len(obs_spec))
    for i, spec in enumerate(obs_spec):
        idx = None
        for k, nm in enumerate(names_list):
            if nm == spec['name']:
                idx = k
                break
        if idx is None:
            raise ValueError("Element '%s' not found in the LQN." % spec['name'])
        metric = str(spec['metric']).lower()
        if metric == 'respt':
            z[i] = metrics['RespT'][idx]
        elif metric == 'util':
            z[i] = metrics['Util'][idx]
        elif metric == 'tput':
            z[i] = metrics['Tput'][idx]
        elif metric == 'qlen':
            z[i] = metrics['QLen'][idx]
        else:
            raise ValueError("Unknown metric '%s'." % spec['metric'])
    return z


def infer_lqn_jacobian(hfun, a, fd_step=1e-3, fd_floor=1e-6):
    """Forward finite-difference sensitivity matrix H = dh/da and h0 = hfun(a).

    Column i is (hfun(a + d_i) - h0) / d_i with d_i = fd_step*max(|a_i|, fd_floor).
    This is the approximate sensitivity matrix H_k used in the EKF update.
    """
    a = np.asarray(a, dtype=float).ravel()
    h0 = np.asarray(hfun(a), dtype=float).ravel()
    H = np.zeros((h0.size, a.size))
    for i in range(a.size):
        d = fd_step * max(abs(a[i]), fd_floor)
        ap = a.copy()
        ap[i] += d
        hi = np.asarray(hfun(ap), dtype=float).ravel()
        H[:, i] = (hi - h0) / d
    return H, h0


def infer_lqn_ekf(hfun, a0, P0, Z, Q, R, options=None):
    """Extended Kalman Filter tracking a hidden parameter vector across Z.

    hfun maps a parameter vector to a predicted observation vector z = h(a).
    a0 (np,), P0 (np,np), Z (no,nsteps), Q (np,np), R (no,no). Returns
    (ahat (np,nsteps), info) where info has keys P, Phist, e, zpred, Er and Ea
    (parameter tracking RMS vs options['aTrue'] when supplied, else None).
    """
    if options is None:
        options = {}
    fd_step = options.get('fdStep', 1e-3)
    fd_floor = options.get('fdFloor', 1e-6)
    clamp = options.get('clampPositive', True)
    a_true = options.get('aTrue', None)
    verbose = options.get('verbose', False)

    a0 = np.asarray(a0, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    no, nsteps = Z.shape
    npar = a0.size

    ahat = np.zeros((npar, nsteps))
    e_hist = np.zeros((no, nsteps))
    zpred_hist = np.zeros((no, nsteps))
    Phist = []

    a = a0.copy()
    P = np.asarray(P0, dtype=float).copy()
    I_np = np.eye(npar)
    for k in range(nsteps):
        # (1) predict: zero-mean drift; project covariance (eq 5)
        a_pred = a
        Ppred = P + Q
        # (2,4) predicted measurement and sensitivity matrix
        H, zpred = infer_lqn_jacobian(hfun, a_pred, fd_step, fd_floor)
        # (3) prediction error
        e = Z[:, k] - zpred
        # (6) Kalman gain (suboptimal because h is nonlinear)
        S = H @ Ppred @ H.T + R
        K = Ppred @ H.T @ np.linalg.inv(S)
        # (4-update) improved estimate
        a = a_pred + K @ e
        if clamp:
            a = np.maximum(a, fd_floor)
        # (7) covariance update; symmetrize for stability
        P = (I_np - K @ H) @ Ppred
        P = 0.5 * (P + P.T)

        ahat[:, k] = a
        e_hist[:, k] = e
        zpred_hist[:, k] = zpred
        Phist.append(P.copy())
        if verbose:
            print('[infer_lqn_ekf] step %d/%d  ||e||=%.4g' % (k + 1, nsteps, np.linalg.norm(e)))

    info = {'P': P, 'Phist': Phist, 'e': e_hist, 'zpred': zpred_hist}
    info['Er'] = float(np.sqrt(np.mean(e_hist ** 2)))
    if a_true is not None:
        a_true = np.asarray(a_true, dtype=float).ravel()
        D = ahat - a_true.reshape(-1, 1)
        info['Ea'] = float(np.sqrt(np.mean(D ** 2)))
    else:
        info['Ea'] = None
    return ahat, info


def infer_lqn(model, param_spec, obs_spec, Z, options=None):
    """Identify hidden LQN parameters from measured performance data.

    Estimates the parameters named in param_spec (activity host demands and/or
    task think times) of the LayeredNetwork model from the measurement sequence
    Z (no x nsteps) using an EKF over the observation model defined by obs_spec.
    A single measurement column with options['QFac'] = 0 reduces to one-shot
    least-squares calibration.

    options keys (defaults): solver (callable model->solver, default native
    SolverLN), QFac (0.1), RFac (0.2), cvA (1), gammaT ([]; else T/Tstar; else 1),
    T, Tstar, Q, R, P0, a0, aTrue, fdStep (1e-3), fdFloor (1e-6),
    clampPositive (True), verbose (False).

    Returns (model, info): model has the final estimate applied; info is the EKF
    result dict augmented with ahat, a0, Q, R, P0.
    """
    if options is None:
        options = {}
    eps = np.finfo(float).eps

    no = len(obs_spec)
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z.reshape(-1, 1)
    if Z.shape[0] != no:
        raise ValueError('Row count of Z must equal len(obs_spec).')

    solver_ctor = options.get('solver', None)

    a0 = options.get('a0', None)
    if a0 is None:
        a0 = infer_lqn_getparams(model, param_spec)
    a0 = np.asarray(a0, dtype=float).ravel()

    QFac = options.get('QFac', 0.1)
    RFac = options.get('RFac', 0.2)
    cvA = options.get('cvA', 1.0)
    gammaT = options.get('gammaT', None)
    if gammaT is None:
        T = options.get('T', None)
        Tstar = options.get('Tstar', None)
        gammaT = (T / Tstar) if (T is not None and Tstar is not None) else 1.0

    Q = options.get('Q', None)
    if Q is None:
        qd = (QFac * np.abs(a0) * cvA) ** 2                 # eq 9a
        Q = np.diag(np.maximum(qd, eps))
    R = options.get('R', None)
    if R is None:
        zbar = np.mean(Z, axis=1)
        rd = ((RFac * np.abs(zbar)) / 1.96) ** 2 / gammaT   # eq 9b
        R = np.diag(np.maximum(rd, eps))
    P0 = options.get('P0', None)
    if P0 is None:
        P0 = np.diag(np.maximum((0.5 * np.abs(a0)) ** 2, eps))

    def hfun(a):
        infer_lqn_setparams(model, param_spec, a)
        if solver_ctor is not None:
            solver = solver_ctor(model)
        else:
            from line_solver import SolverLN
            solver = SolverLN(model, verbose=False)
        QN, UN, RN, TN, _AN, _WN = solver.get_ensemble_avg()
        names = solver.lqn.names
        metrics = {'QLen': QN, 'Util': UN, 'RespT': RN, 'Tput': TN}
        return infer_lqn_getobs(names, metrics, obs_spec)

    ekfopts = {
        'fdStep': options.get('fdStep', 1e-3),
        'fdFloor': options.get('fdFloor', 1e-6),
        'clampPositive': options.get('clampPositive', True),
        'aTrue': options.get('aTrue', None),
        'verbose': options.get('verbose', False),
    }
    ahat, info = infer_lqn_ekf(hfun, a0, P0, Z, Q, R, ekfopts)
    info['ahat'] = ahat
    info['a0'] = a0
    info['Q'] = Q
    info['R'] = R
    info['P0'] = P0

    infer_lqn_setparams(model, param_spec, ahat[:, -1])
    return model, info
