import numpy as np
from line_solver import SolverMVA


def infer_generate_qlen_traces(model, stride):
    """Generate synthetic queue length trace data.

    Uses MVA steady-state queue lengths with added noise to simulate
    transient traces. Python LINE SSA does not support transient traces,
    so this generates synthetic data for testing purposes.

    Args:
        model: A LINE Network model.
        stride: Stride for downsampling transient results.

    Returns:
        Tuple (timeIntervals, queueLengthMatrix) where:
            timeIntervals: Downsampled time points array.
            queueLengthMatrix: Nested list [station][class] of queue length metric arrays.

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    solver = SolverMVA(model)
    avgQLen = solver.getAvgQLen()

    station_count = model.getNumberOfNodes()
    job_count = len(model.classes)

    n_points = 1000
    time_intervals = np.linspace(0, 5, n_points)
    time_intervals = time_intervals[::stride]
    n_out = len(time_intervals)

    avgQLengths = [[None for _ in range(job_count)] for _ in range(station_count)]

    for i in range(station_count):
        for c in range(job_count):
            if avgQLen.ndim == 2:
                ql_mean = avgQLen[i, c]
            else:
                ql_mean = avgQLen[i]
            # Generate noisy trace around steady-state
            noise = np.random.randn(n_out) * 0.01 * max(ql_mean, 1.0)
            avgQLengths[i][c] = np.maximum(0, ql_mean + noise)

    return time_intervals, avgQLengths
