"""
Polling System Analysis Algorithms.

Native Python implementations for analyzing polling/vacation queue systems
with various disciplines including exhaustive, gated, and 1-limited.

References:
    Takagi, H., "Queuing analysis of polling models: An update."
    ACM Computing Surveys, Vol. 20, No. 1, March 1988.
"""

import numpy as np
from scipy import linalg
from typing import List, Tuple, Union

# Import MAP functions from our native MAM module
from ..mam import map_lambda, map_mean, map_moment, map_var


def polling_qsys_exhaustive(
    arvMAPs: List[Tuple[np.ndarray, np.ndarray]],
    svcMAPs: List[Tuple[np.ndarray, np.ndarray]],
    switchMAPs: List[Tuple[np.ndarray, np.ndarray]]
) -> np.ndarray:
    """
    Compute exact mean waiting times for exhaustive polling system.

    In exhaustive polling, the server continues to serve a queue until it
    becomes empty before moving to the next queue.

    Based on Takagi, ACM Computing Surveys, Vol. 20, No. 1, 1988, eq (15).

    Args:
        arvMAPs: List of arrival process MAPs, each as (D0, D1) tuple.
        svcMAPs: List of service process MAPs, each as (D0, D1) tuple.
        switchMAPs: List of switching time MAPs, each as (D0, D1) tuple.

    Returns:
        Array of mean waiting times for each queue.
    """
    n = len(arvMAPs)

    # Extract MAP parameters
    lambda_arr = np.zeros(n)
    b = np.zeros(n)
    b2 = np.zeros(n)
    rho1 = np.zeros(n)
    r1 = np.zeros(n)
    delta2 = np.zeros(n)
    rho = 0.0
    r = 0.0

    for i in range(n):
        D0, D1 = arvMAPs[i]
        lambda_arr[i] = map_lambda(D0, D1)

        D0, D1 = svcMAPs[i]
        b[i] = map_mean(D0, D1)
        b2[i] = map_moment(D0, D1, 2)

        rho1[i] = lambda_arr[i] * b[i]
        rho += rho1[i]

        D0, D1 = switchMAPs[i]
        r1[i] = map_mean(D0, D1)
        r += r1[i]
        delta2[i] = map_var(D0, D1)

    # Build system of equations
    lst1 = []
    lst2 = []

    for i in range(1, n + 1):  # 1-indexed as in Kotlin
        for j in range(1, n + 1):
            if i > j:
                t1 = np.zeros(n * n)
                for m in range(i + 1, n + 1):
                    t1[(j - 1) * n + m - 1] -= 1
                for m in range(1, j):
                    t1[(j - 1) * n + m - 1] -= 1
                for m in range(j, i):
                    t1[(m - 1) * n + j - 1] -= 1
                t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1]
                lst1.append(t1)
                lst2.append(0.0)
            elif j > i:
                t1 = np.zeros(n * n)
                for m in range(i + 1, j):
                    t1[(j - 1) * n + m - 1] -= 1
                for m in range(j, n + 1):
                    t1[(m - 1) * n + j - 1] -= 1
                for m in range(1, i):
                    t1[(m - 1) * n + j - 1] -= 1
                t1[(i - 1) * n + j - 1] += (1 - rho1[i - 1]) / rho1[i - 1]
                lst1.append(t1)
                lst2.append(0.0)
            else:  # i == j
                t1 = np.zeros(n * n)
                t1[(i - 1) * n + i - 1] += 1
                for m in range(1, n + 1):
                    if i != m:
                        t1[(i - 1) * n + m - 1] -= rho1[i - 1] / (1 - rho1[i - 1])
                lst1.append(t1)

                temp = 0.0
                idx = i - 2 if i > 1 else n - 1
                temp += delta2[idx] / (1 - rho1[i - 1])**2
                temp += (lambda_arr[i - 1] * b2[i - 1] * r * (1 - rho1[i - 1]) /
                         ((1 - rho) * (1 - rho1[i - 1])**3))
                lst2.append(temp)

    # Solve linear system
    size = len(lst1)
    A = np.array(lst1)  # Shape: (size, n*n)
    B = np.array(lst2)  # Shape: (size,)

    try:
        x = linalg.lstsq(A, B)[0]
    except:
        x = np.zeros(n * n)

    # Compute waiting times
    W = np.zeros(n)
    for i in range(1, n + 1):
        temp = lambda_arr[i - 1] * b2[i - 1] / (2 * (1 - rho1[i - 1]))
        temp += r * (1 - rho1[i - 1]) / (2 * (1 - rho))

        sum_val = 0.0
        for j in range(1, n + 1):
            if i != j:
                sum_val += x[(i - 1) * n + j - 1]

        sum_val *= (1 - rho1[i - 1]) / rho1[i - 1]
        idx = i - 2 if i > 1 else n - 1
        sum_val += delta2[idx]
        sum_val /= r * (1 - rho1[i - 1]) * 2 / (1 - rho)
        temp += sum_val

        W[i - 1] = temp

    return W


def polling_qsys_gated(
    arvMAPs: List[Tuple[np.ndarray, np.ndarray]],
    svcMAPs: List[Tuple[np.ndarray, np.ndarray]],
    switchMAPs: List[Tuple[np.ndarray, np.ndarray]]
) -> np.ndarray:
    """
    Compute exact mean waiting times for gated polling system.

    In gated polling, the server serves all customers present at the
    beginning of a visit period.

    Based on Takagi, ACM Computing Surveys, Vol. 20, No. 1, 1988, eq (20).

    Args:
        arvMAPs: List of arrival process MAPs, each as (D0, D1) tuple.
        svcMAPs: List of service process MAPs, each as (D0, D1) tuple.
        switchMAPs: List of switching time MAPs, each as (D0, D1) tuple.

    Returns:
        Array of mean waiting times for each queue.
    """
    n = len(arvMAPs)

    # Extract MAP parameters
    lambda_arr = np.zeros(n)
    b = np.zeros(n)
    b2 = np.zeros(n)
    rho1 = np.zeros(n)
    r1 = np.zeros(n)
    delta2 = np.zeros(n)
    rho = 0.0
    r = 0.0

    for i in range(n):
        D0, D1 = arvMAPs[i]
        lambda_arr[i] = map_lambda(D0, D1)

        D0, D1 = svcMAPs[i]
        b[i] = map_mean(D0, D1)
        b2[i] = map_moment(D0, D1, 2)

        rho1[i] = lambda_arr[i] * b[i]
        rho += rho1[i]

        D0, D1 = switchMAPs[i]
        r1[i] = map_mean(D0, D1)
        r += r1[i]
        delta2[i] = map_var(D0, D1)

    # Build system of equations
    lst1 = []
    lst2 = []

    for i in range(n):  # 0-indexed
        for j in range(n):
            if i > j:
                t1 = np.zeros(n * n)
                for m in range(i, n):
                    t1[j * n + m] = -1.0
                for m in range(0, j):
                    t1[j * n + m] = -1.0
                for m in range(j, i):
                    t1[m * n + j] = -1.0
                t1[i * n + j] = 1 / rho1[i]
                lst1.append(t1)
                lst2.append(0.0)
            elif j > i:
                t1 = np.zeros(n * n)
                for m in range(i, j):
                    t1[j * n + m] = -1.0
                for m in range(j, n):
                    t1[m * n + j] = -1.0
                for m in range(0, i):
                    t1[m * n + j] = -1.0
                t1[i * n + j] = 1 / rho1[i]
                lst1.append(t1)
                lst2.append(0.0)
            else:  # i == j
                t1 = np.zeros(n * n)
                t1[i * n + i] += 1
                for m in range(n):
                    if i != m:
                        t1[i * n + m] = -rho1[i]
                for m in range(n):
                    t1[m * n + i] -= rho1[i] * rho1[i]
                lst1.append(t1)
                temp = delta2[i] + lambda_arr[i] * b2[i] * r / (1 - rho)
                lst2.append(temp)

    # Solve linear system
    n1 = len(lst1)
    A = np.array(lst1)
    B = np.array(lst2)

    try:
        x = linalg.lstsq(A, B)[0]
    except:
        x = np.zeros(n * n)

    # Compute waiting times
    W = np.zeros(n)
    for i in range(n):
        temp = (1 + rho1[i]) * r / (2 * (1 - rho))

        sum_val = 0.0
        for j in range(n):
            if i != j:
                sum_val += x[i * n + j]
        sum_val *= (1 / rho1[i])

        for j in range(n):
            sum_val += x[j * n + i]

        temp += (1 - rho) * (1 + rho1[i]) * sum_val / (2 * r)
        W[i] = temp

    return W


def polling_qsys_1limited(
    arvMAPs: List[Tuple[np.ndarray, np.ndarray]],
    svcMAPs: List[Tuple[np.ndarray, np.ndarray]],
    switchMAPs: List[Tuple[np.ndarray, np.ndarray]]
) -> np.ndarray:
    """
    Compute exact mean waiting times for 1-limited polling system.

    In 1-limited polling, the server serves at most one customer from each
    queue before moving to the next queue.

    Based on Takagi, ACM Computing Surveys, Vol. 20, No. 1, 1988, eq (20).

    Args:
        arvMAPs: List of arrival process MAPs, each as (D0, D1) tuple.
        svcMAPs: List of service process MAPs, each as (D0, D1) tuple.
        switchMAPs: List of switching time MAPs, each as (D0, D1) tuple.

    Returns:
        Array of mean waiting times for each queue.
    """
    n = len(arvMAPs)

    # Extract MAP parameters
    lambda_arr = np.zeros(n)
    b = np.zeros(n)
    b2 = np.zeros(n)
    rho1 = np.zeros(n)
    r1 = np.zeros(n)
    delta2 = np.zeros(n)
    rho = 0.0
    r = 0.0
    d = 0.0

    for i in range(n):
        D0, D1 = arvMAPs[i]
        lambda_arr[i] = map_lambda(D0, D1)

        D0, D1 = svcMAPs[i]
        b[i] = map_mean(D0, D1)
        b2[i] = map_moment(D0, D1, 2)

        rho1[i] = lambda_arr[i] * b[i]
        rho += rho1[i]

        D0, D1 = switchMAPs[i]
        r1[i] = map_mean(D0, D1)
        r += r1[i]
        delta2[i] = map_var(D0, D1)
        d += delta2[i]

    # Compute waiting times using closed-form formula
    W = np.zeros(n)

    sum_of_squares = np.sum(rho1**2)
    sum_product = np.sum(lambda_arr * b2)

    for i in range(n):
        W[i] = (1 - rho + rho1[i]) / (1 - rho - lambda_arr[i] * r)
        W[i] *= (1 - rho) / ((1 - rho) * rho + sum_of_squares)
        W[i] *= (rho / (2 * (1 - rho)) * sum_product +
                 rho * d / 2 / r +
                 r / (2 * (1 - rho)) * rho1[i] * (1 + rho1[i]))

    return W


__all__ = [
    'polling_qsys_exhaustive',
    'polling_qsys_gated',
    'polling_qsys_1limited',
]
