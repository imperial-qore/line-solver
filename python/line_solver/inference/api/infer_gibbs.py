import numpy as np

from line_solver.api.pfqn.mva import pfqn_bs


def infer_gibbs(data, nb_cores, tol=1e-3):
    """Gibbs sampling demand estimation from trace data.

    Args:
        data: list of lists in standard format.
              data[2][k] = arrival times (ms), data[3][k] = response times,
              data[5][k] = throughput
        nb_cores: number of server cores
        tol: convergence tolerance

    Returns:
        demand: 1-D array of estimated demands

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    data_needed = 200000
    likelihood_sample = 5000
    nb_samples = 2000

    nb_classes = len(data[2])
    nb_nodes = 2
    nb_jobs = np.zeros(nb_classes, dtype=int)

    prob, nb_jobs, N0 = _analyse_data(data, nb_jobs, nb_classes, nb_nodes, data_needed)

    used_cores = 0.0
    for k in range(prob.shape[0]):
        queue_jobs = np.sum(prob[k, nb_classes:2 * nb_classes])
        if queue_jobs > nb_cores:
            used_cores += nb_cores * prob[k, -1]
        else:
            used_cores += queue_jobs * prob[k, -1]
    used_cores = used_cores / (1 - prob[-1, -1])

    think_time = np.zeros(nb_classes)
    for k in range(nb_classes):
        tput = np.mean(data[5][k])
        if tput > 0:
            think_time[k] = (nb_jobs[k] - N0[k]) / tput

    range_size = np.ones(nb_classes * (nb_nodes - 1))

    cum_prob = np.cumsum(prob[:, nb_classes * nb_nodes])
    testset = np.zeros((likelihood_sample, nb_classes * nb_nodes))
    for k in range(likelihood_sample):
        uni_value = np.random.rand()
        index = np.searchsorted(cum_prob, uni_value)
        if index >= prob.shape[0]:
            index = prob.shape[0] - 1
        testset[k, :] = prob[index, :nb_classes * nb_nodes]

    # Log factorial table
    LV = np.zeros(int(np.sum(nb_jobs)) + 2)
    for k in range(2, len(LV)):
        LV[k] = LV[k - 1] + np.log(k - 1)

    sum_a = np.sum(LV[(testset.astype(int)).flatten() + 1])

    initial = np.zeros(nb_classes * (nb_nodes - 1))
    log_g_initial = np.sum(nb_jobs * np.log(think_time + 1e-300))
    for k in range(nb_classes):
        log_g_initial -= np.sum(np.log(np.arange(1, nb_jobs[k] + 1)))

    smpl = np.zeros((nb_samples, nb_classes * (nb_nodes - 1)))
    sample_index = 0

    for k in range(nb_samples // 50):
        for s in range(50):
            for h in range(nb_classes * (nb_nodes - 1)):
                if sample_index == 0:
                    theta = np.concatenate([smpl[sample_index, :h], initial[h:]])
                else:
                    theta = np.concatenate([smpl[sample_index, :h], smpl[sample_index - 1, h:]])

                smpl[sample_index, h], log_g_initial, rs_dim = _gibbs_sampler_simple(
                    think_time, theta, testset, h, nb_nodes, nb_classes,
                    nb_jobs, log_g_initial, tol, range_size[h], LV, sum_a)
                range_size[h] = rs_dim * 2

            sample_index += 1

        if k == 1:
            demand_old = np.mean(smpl[50:sample_index, :], axis=0)
        elif k > 1:
            demand_now = np.mean(smpl[k * 50:(k + 1) * 50, :], axis=0)
            demand_now = demand_now / (k + 1) + demand_old / (k + 1) * k
            if np.mean(np.abs((demand_now - demand_old) / (demand_old + 1e-300))) < tol:
                N_cut = sample_index // 2
                demand = np.mean(smpl[N_cut:sample_index, :] * used_cores, axis=0)
                return demand
            else:
                demand_old = demand_now

    N_cut = sample_index // 2
    demand = np.mean(smpl[N_cut:sample_index, :] * used_cores, axis=0)
    return demand


def _analyse_data(data, nb_jobs, nb_classes, nb_nodes, data_needed):
    K = nb_classes
    N = nb_jobs.copy()
    N0 = np.zeros(K)

    temp_ts = []
    temp_class = []
    temp_logger = []
    for i in range(K):
        arr = np.asarray(data[2][i], dtype=float).flatten()
        rt = np.asarray(data[3][i], dtype=float).flatten()
        temp_length = len(arr)
        temp_ts.extend(arr.tolist())
        temp_ts.extend((arr + rt * 1000).tolist())
        temp_class.extend([i + 1] * temp_length * 2)
        temp_logger.extend([1] * temp_length)
        temp_logger.extend([2] * temp_length)

    temp_ts = np.array(temp_ts)
    temp_class = np.array(temp_class, dtype=int)
    temp_logger = np.array(temp_logger, dtype=int)

    order = np.argsort(temp_ts, kind='stable')
    ts = temp_ts[order]
    class_id = temp_class[order]
    logger_id = temp_logger[order]

    burnin = max(0, len(ts) - data_needed)
    if data_needed == 0:
        burnin = 0

    total_length = len(ts)
    count = np.zeros((total_length, K, nb_nodes))
    count[0, :, 0] = N

    for i in range(total_length - 1):
        count[i + 1, :, :] = count[i, :, :]
        c = class_id[i] - 1
        l = logger_id[i] - 1
        count[i + 1, c, l] -= 1
        if l == nb_nodes - 1:
            count[i + 1, c, 0] += 1
        else:
            count[i + 1, c, l + 1] += 1

    if np.sum(N) == 0:
        for i in range(K):
            N[i] = int(np.max(count[:, i, :]))

    for i in range(total_length):
        for j in range(K):
            count[i, j, 0] += N[j]

    count_2d = count.reshape(total_length, K * nb_nodes)
    time_interval = np.zeros(total_length)
    time_interval[1:] = np.diff(ts)

    count_with_time = np.column_stack([count_2d, time_interval])
    count_with_time = count_with_time[burnin:]

    state_data = count_with_time[:, :-1]
    time_data = count_with_time[:, -1]

    unique_rows, inverse = np.unique(state_data, axis=0, return_inverse=True)
    time_col = np.zeros(len(unique_rows))
    for i in range(len(unique_rows)):
        time_col[i] = np.sum(time_data[inverse == i])

    obs_length = ts[-1] - ts[burnin]
    prob = np.column_stack([unique_rows, time_col / obs_length])

    for i in range(K):
        N0[i] = np.sum(prob[:, -1] * prob[:, K + i])

    return prob, N, N0


def _gibbs_sampler_simple(think_time, theta, testset, index, nb_nodes,
                          nb_classes, nb_jobs, log_g_initial, interval,
                          range_size, LV, sum_a):
    steps = np.arange(0, range_size + interval, interval)
    N_steps = len(steps)

    x = np.zeros((nb_nodes, nb_classes))
    x[0, :] = think_time
    for i in range(nb_nodes - 1):
        x[i + 1, :] = theta[i * nb_classes:(i + 1) * nb_classes]

    index_i = index // nb_classes + 1  # 1-based node index (queue)
    index_j = index % nb_classes  # 0-based class index

    log_g = np.zeros(N_steps)

    _, QN, _, _, _ = pfqn_bs(x[1:, :], nb_jobs, x[0, :])

    index_previous = np.searchsorted(steps, theta[index])
    if index_previous >= N_steps:
        index_previous = N_steps - 1
    log_g[index_previous] = log_g_initial

    for i in range(index_previous - 1, -1, -1):
        x[index_i, index_j] = steps[i + 1]
        _, QN, _, _, _ = pfqn_bs(x[1:, :], nb_jobs, x[0, :], interval, 1000, QN)
        val = 1 + QN[index_i - 1, index_j] / (steps[i + 1] + 1e-300) * (-interval)
        if val > 0:
            log_g[i] = log_g[i + 1] + np.log(val)
        else:
            log_g[i] = log_g[i + 1]

    x[index_i, index_j] = theta[index]
    _, QN, _, _, _ = pfqn_bs(x[1:, :], nb_jobs, x[0, :])

    for i in range(index_previous + 1, N_steps):
        x[index_i, index_j] = steps[i - 1]
        _, QN, _, _, _ = pfqn_bs(x[1:, :], nb_jobs, x[0, :], interval, 1000, QN)
        val = 1 + QN[index_i - 1, index_j] / (steps[i - 1] + 1e-300) * interval
        if val > 0:
            log_g[i] = log_g[i - 1] + np.log(val)
        else:
            log_g[i] = log_g[i - 1]

    log_prob = np.zeros(N_steps)
    for i in range(N_steps):
        theta[index] = steps[i]
        log_prob[i] = np.sum(testset[:, index + nb_classes]) * np.log(steps[i] + 1e-300) \
                       - log_g[i] * testset.shape[0]

    log_prob -= np.max(log_prob)
    prob = np.exp(log_prob)
    prob /= np.sum(prob)

    cum_prob = np.cumsum(prob)
    rand_var = np.random.rand()
    index_prob = np.searchsorted(cum_prob, rand_var)

    range_size_dim = np.searchsorted(cum_prob, 1 - 1e-10)
    if range_size_dim < N_steps:
        range_size_dim = steps[range_size_dim] * 2
    else:
        range_size_dim = steps[-1] * 2

    if index_prob >= N_steps:
        value = theta[index]
        log_g_current = log_g_initial
    else:
        value = steps[index_prob]
        log_g_current = log_g[index_prob]

    return value, log_g_current, range_size_dim
