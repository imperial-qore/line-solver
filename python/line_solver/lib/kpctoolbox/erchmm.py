"""
Extended Renewal Continuous-time Hidden Markov Model (ER-CHMM) functions.

Ported from MATLAB: matlab/lib/kpctoolbox/erchmm/
JAR reference: jar/src/main/kotlin/jline/lib/kpctoolbox/erchmm/ERCHMM.kt

The EM algorithm is based on:
- Okamura et al. (2008). An EM algorithm for a Superposition of Markovian
  Arrival Processes.
- Horvath & Okamura (2013). A Fast EM algorithm for Fitting Marked Markov
  Arrival Processes with a new Special Structure.
- Horvath et al. (2018). Parallel Algorithms for Fitting Markov Arrival
  Processes.

This algorithm was inspired and adapted from the BUTools implementation:
http://webspn.hit.bme.hu/~telek/tools/butools/doc/MAPFromTrace.html
"""

import numpy as np
from math import factorial, log, log2
import time
import sys


def _all_erlang(branches, sum_erlangs):
    """
    Find all unique combinations of Erlang branch orders that sum to
    sum_erlangs using the given number of branches.

    Each branch must have order >= 1. Combinations are stored sorted
    to avoid duplicates.

    Parameters
    ----------
    branches : int
        Number of Erlang branches.
    sum_erlangs : int
        Total sum of orders across all branches.

    Returns
    -------
    list of list of int
        List of unique sorted combinations.
    """
    if branches == 1:
        return [[sum_erlangs]]

    erlangs = []
    for k1 in range(1, sum_erlangs - branches + 2):
        sub_combinations = _all_erlang(branches - 1, sum_erlangs - k1)
        for sub in sub_combinations:
            sorted_erlang = sorted(sub + [k1])
            # Check if this combination already exists
            erlang_found = False
            for existing in erlangs:
                if existing == sorted_erlang:
                    erlang_found = True
                    break
            if not erlang_found:
                erlangs.append(sorted_erlang)
    return erlangs


def _generate_D0_from_erlangs(lambda_vals, orders):
    """
    Generate the D0 matrix from Erlang representation.

    Each Erlang branch i has order orders[i] and rate lambda_vals[i].
    The resulting matrix is block-diagonal with each block being an
    Erlang sub-generator of the corresponding order and rate.

    Parameters
    ----------
    lambda_vals : array_like
        Rate parameters for each Erlang branch, shape (M,).
    orders : array_like of int
        Orders of each Erlang branch, shape (M,).

    Returns
    -------
    numpy.ndarray
        D0 matrix of shape (sum(orders), sum(orders)).
    """
    n = int(np.sum(orders))
    D0 = np.zeros((n, n))
    init_x = 0
    for i in range(len(lambda_vals)):
        order_i = int(orders[i])
        lam = lambda_vals[i]
        # Build the Erlang sub-generator block:
        # diagonal = -lambda, superdiagonal = +lambda
        block = lam * (np.diag(np.ones(order_i - 1), 1) - np.diag(np.ones(order_i)))
        D0[init_x:init_x + order_i, init_x:init_x + order_i] = block
        init_x += order_i
    return D0


def erchmm_emfit(trace, orders, iter_max=300, iter_tol=1e-7, verbose=True):
    """
    Fit an Extended Renewal Continuous-time Hidden Markov Model (ER-CHMM)
    to a trace using the EM algorithm.

    When `orders` is a scalar (single integer), the algorithm searches over
    all unique combinations of Erlang branch orders summing to that value
    and returns the MAP with the best log-likelihood.

    When `orders` is a list/array of integers, the algorithm directly fits
    an ER-CHMM with those branch orders.

    Parameters
    ----------
    trace : array_like
        Array of inter-arrival times (positive reals), shape (K,).
    orders : int or array_like of int
        If a single integer: total sum of Erlang orders; all combinations
        of branches are tried.
        If a list/array: the Erlang order for each branch.
    iter_max : int, optional
        Maximum number of EM iterations (default 300).
    iter_tol : float, optional
        Convergence tolerance on relative change of log-likelihood
        (default 1e-7).
    verbose : bool, optional
        Whether to print progress messages (default True).

    Returns
    -------
    MAP : list of numpy.ndarray
        Two-element list [D0, D1] representing the fitted MAP.
    logL : float
        Log-likelihood of the fitted model (per observation).
    """
    trace = np.asarray(trace, dtype=float).ravel()

    # Handle scalar orders: try all combinations
    if np.isscalar(orders) or (hasattr(orders, '__len__') and len(orders) == 1):
        total_order = int(orders) if np.isscalar(orders) else int(orders[0])

        X_chosen = None
        Y_chosen = None
        orders_best = None
        log_li_best = -np.inf

        for ord_num_1 in range(2, total_order + 1):
            all_orders = _all_erlang(ord_num_1, total_order)
            for ord_num_2 in range(len(all_orders)):
                if verbose:
                    orders_str = ','.join(str(o) for o in all_orders[ord_num_2])
                    print(f'Calculating with orders {orders_str}...')
                    sys.stdout.flush()

                # Recursive call with specific orders
                ord_MAP, l = erchmm_emfit(
                    trace, all_orders[ord_num_2], iter_max, iter_tol, verbose
                )
                ord_X = ord_MAP[0]
                ord_Y = ord_MAP[1]
                if l > log_li_best:
                    X_chosen = ord_X
                    Y_chosen = ord_Y
                    log_li_best = l
                    orders_best = all_orders[ord_num_2]

        # Finalize
        D0 = X_chosen
        D1 = Y_chosen
        logL = log_li_best

        if verbose:
            orders_str = ','.join(str(o) for o in orders_best) if orders_best is not None else ''
            print(f'Best solution: log-likelihood={log_li_best}, orders={orders_str}')
            sys.stdout.flush()

        return [D0, D1], logL

    # --- Main EM algorithm for a given orders vector ---

    orders = np.asarray(orders, dtype=int).ravel()
    M = len(orders)
    K = len(trace)

    # Initialize pi and lambda parameters such that the mean is matched
    pi_v = np.ones(M) / M
    # MATLAB: lambda = diag(orders) * (1:M)'
    lambda_vals = orders.astype(float) * np.arange(1, M + 1, dtype=float)
    trace_mean = np.sum(trace) / K
    pi_mean = np.sum(pi_v / np.arange(1, M + 1, dtype=float))
    # Match the mean
    lambda_vals = lambda_vals * pi_mean / trace_mean

    # Initialize the transition matrix T = ones(M,1) * pi_v
    T = np.outer(np.ones(M), pi_v)

    # Initialize matrices for the EM algorithm
    F = np.zeros((M, K))                # Branch densities
    A_likelihoods = np.zeros((K, M))    # Forward likelihoods
    B_likelihoods = np.zeros((M, K))    # Backward likelihoods
    A_likelihoods_scale = np.zeros(K)   # Forward scale factors
    B_likelihoods_scale = np.zeros(K)   # Backward scale factors

    ologli = 1.0
    logL = 0.0
    steps = 1
    t1 = time.time()

    # Main EM iteration loop
    # MATLAB: abs((1-0)/0) = Inf > tol, so the loop enters on first iteration.
    # In Python we must guard against division by zero.
    while (logL == 0.0 or abs((ologli - logL) / logL) > iter_tol) and steps < iter_max:
        ologli = logL

        # --- E-step ---

        # Compute branch densities
        for i in range(M):
            # F(i,:) = ((lambda(i)*trace).^(orders(i)-1) / factorial(orders(i)-1)
            #           * lambda(i)) .* exp(-lambda(i)*trace)
            lt = lambda_vals[i] * trace
            order_i = int(orders[i])
            F[i, :] = (
                (lt ** (order_i - 1) / factorial(order_i - 1) * lambda_vals[i])
                * np.exp(-lambda_vals[i] * trace)
            )

        # Compute forward likelihood vectors
        prev_pi = pi_v.copy()
        scaled_prev = 0.0
        for k in range(K):
            # prev_pi = prev_pi * diag(F(:,k)) * T
            prev_pi = (prev_pi * F[:, k]) @ T
            s = np.sum(prev_pi)
            if s > 0:
                scale = log2(s)
            else:
                scale = 0.0
            prev_pi = prev_pi * 2.0 ** (-scale)
            A_likelihoods_scale[k] = scaled_prev + scale
            A_likelihoods[k, :] = prev_pi.copy()
            scaled_prev = A_likelihoods_scale[k]

        # a_forward_likelihoods = [pi_v; A_likelihoods(1:end-1,:)]
        a_forward_likelihoods = np.zeros((K, M))
        a_forward_likelihoods[0, :] = pi_v
        a_forward_likelihoods[1:, :] = A_likelihoods[:-1, :]
        # A_scaled_v = [0, A_likelihoods_scale(1:end-1)]
        A_scaled_v = np.zeros(K)
        A_scaled_v[1:] = A_likelihoods_scale[:-1]

        # Compute backward likelihood vectors
        next_b = np.ones(M)
        scaled_prev = 0.0
        for k in range(K - 1, -1, -1):
            # next_b = diag(F(:,k)) * T * next_b
            # next_b_new[i] = F[i,k] * sum_j(T[i,j] * next_b[j])
            next_b = F[:, k] * (T @ next_b)
            s = np.sum(next_b)
            if s > 0:
                scale = log2(s)
            else:
                scale = 0.0
            next_b = next_b * 2.0 ** (-scale)
            B_likelihoods_scale[k] = scaled_prev + scale
            B_likelihoods[:, k] = next_b.copy()
            scaled_prev = B_likelihoods_scale[k]

        # b_backward_likelihoods = [B_likelihoods(:,2:end), ones(M,1)]
        b_backward_likelihoods = np.zeros((M, K))
        b_backward_likelihoods[:, :-1] = B_likelihoods[:, 1:]
        b_backward_likelihoods[:, -1] = 1.0
        # B_scale_v = [B_likelihoods_scale(2:end), 0]
        B_scale_v = np.zeros(K)
        B_scale_v[:-1] = B_likelihoods_scale[1:]

        # Compute likelihood
        likelihood_value = pi_v @ B_likelihoods[:, 0]
        if likelihood_value > 0:
            logL = (log(likelihood_value) + B_likelihoods_scale[0] * log(2.0)) / K
        else:
            logL = -np.inf
        i_likelihood = 1.0 / likelihood_value if likelihood_value != 0 else 0.0

        # --- M-step ---

        # Calculate likelihoods_multiplied = a_forward_likelihoods .* B_likelihoods'
        # a_forward_likelihoods is (K, M), B_likelihoods is (M, K), so B_likelihoods.T is (K, M)
        likelihoods_multiplied = a_forward_likelihoods * B_likelihoods.T

        # Normalize each row
        summed_lm = np.sum(likelihoods_multiplied, axis=1, keepdims=True)
        summed_lm[summed_lm == 0] = 1.0  # avoid division by zero
        likelihoods_multiplied = likelihoods_multiplied / summed_lm

        # Compute numerator and denominator for lambda estimations
        numerator_estimation = np.sum(likelihoods_multiplied, axis=0)  # shape (M,)
        denominator_estimation = trace @ likelihoods_multiplied         # shape (M,)

        # Update pi and lambda
        pi_v = numerator_estimation / K
        # Guard against zero denominator
        safe_denom = denominator_estimation.copy()
        safe_denom[safe_denom == 0] = 1.0
        lambda_vals = (orders * numerator_estimation / safe_denom).astype(float)

        # Compute multiplication of forward likelihood and densities
        # dens_mult_a_likelihood = a_forward_likelihoods .* F'
        dens_mult_a_likelihood = a_forward_likelihoods * F.T  # (K, M)

        # Scale factors: summed_lm = i_likelihood * 2.^(A_scaled_v + B_scale_v - B_likelihoods_scale(1))
        summed_lm_vec = i_likelihood * 2.0 ** (A_scaled_v + B_scale_v - B_likelihoods_scale[0])
        # Multiply each row of dens_mult_a_likelihood by the scale factor
        dens_mult_a_likelihood = dens_mult_a_likelihood * summed_lm_vec[:, np.newaxis]

        # Update T: T = (dens_mult_a_likelihood' * b_backward_likelihoods') .* T
        # dens_mult_a_likelihood is (K, M), dens_mult_a_likelihood.T is (M, K)
        # b_backward_likelihoods is (M, K), b_backward_likelihoods.T is (K, M)
        T = (dens_mult_a_likelihood.T @ b_backward_likelihoods.T) * T

        # Normalize T rows
        row_sums = np.sum(T, axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1.0  # avoid division by zero
        T = T / row_sums

        steps += 1

        # Print progress report
        if verbose and (time.time() - t1) > 2:
            print(f'Num of iterations: {steps}, log-likelihood: {logL}')
            sys.stdout.flush()
            t1 = time.time()

    # Show final progress
    if verbose:
        print(f'Num of iterations: {steps}, log-likelihood: {logL}')
        orders_str = ','.join(str(o) for o in orders)
        print(f'EM algorithm terminated. (orders={orders_str})')
        sys.stdout.flush()

    # --- Finalize D0 and D1 matrices ---

    # Generate D0 from lambda and orders
    D0 = _generate_D0_from_erlangs(lambda_vals, orders)

    # Generate D1
    n = int(np.sum(orders))
    D1 = np.zeros((n, n))
    # indicesTo = [1, cumsum(orders(1:end-1))+1] (MATLAB 1-indexed)
    # indicesFrom = cumsum(orders) (MATLAB 1-indexed)
    # Python 0-indexed:
    indices_to = np.zeros(M, dtype=int)
    indices_from = np.zeros(M, dtype=int)
    cumsum_orders = np.cumsum(orders)
    indices_to[0] = 0
    for i in range(1, M):
        indices_to[i] = cumsum_orders[i - 1]
    indices_from = cumsum_orders - 1

    # D1(indicesFrom, indicesTo) = diag(lambda) * T
    for i in range(M):
        for j in range(M):
            D1[indices_from[i], indices_to[j]] = lambda_vals[i] * T[i, j]

    MAP = [D0, D1]
    return MAP, logL
