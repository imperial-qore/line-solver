import numpy as np
from line_solver import SolverSSA


def infer_generate_transient_samples(model, samples, C, stride):
    """Generate transient performance metric samples using SSA simulation.

    Args:
        model: A LINE Network model.
        samples: Number of sample replications for arrival rates.
        C: Number of job classes.
        stride: Stride for downsampling transient results.

    Returns:
        Tuple (arvR, util, avgQLengths, timeIntervals) where:
            arvR: numpy array with replicated average arrival rates.
            util: Nested list [station][class] of utilization metric arrays.
            avgQLengths: Nested list [station][class] of queue length metric arrays.
            timeIntervals: Downsampled time points array.

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    solver = SolverSSA(model)
    solver.getAvg()
    arv_r = solver.getAvgArvR()
    arvR = np.tile(arv_r, (samples, 1))

    Q, U, T = model.getTranHandles()
    QNt, Ut, _ = SolverSSA(model, force=True, timespan=[0, 5]).getTranAvg(Q, U, T)

    station_count = model.getNumberOfNodes()
    job_count = len(model.classes)

    avgQLengths = [[None for _ in range(job_count)] for _ in range(station_count)]
    util = [[None for _ in range(job_count)] for _ in range(station_count)]

    time_intervals = QNt[0][0].t
    time_intervals = time_intervals[::stride]

    for i in range(station_count):
        for c in range(job_count):
            met = QNt[i][c].metric
            avgQLengths[i][c] = met[::stride]
            umet = Ut[i][c].metric
            util[i][c] = umet[::stride]

    return arvR, util, avgQLengths, time_intervals
