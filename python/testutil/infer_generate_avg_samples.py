import numpy as np
from line_solver import SolverSSA


def infer_generate_avg_samples(model, samples, C):
    """Generate average performance metric samples using SSA simulation.

    Args:
        model: A LINE Network model.
        samples: Number of simulation samples to generate.
        C: Number of job classes.

    Returns:
        Tuple (respT, arvR, util) where:
            respT: numpy array of shape (samples, numStations, C) with average response times.
            arvR: numpy array of shape (samples, numStations, C) with average arrival rates.
            util: numpy array of shape (samples, numStations) with total utilizations.

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    num_stations = len(model.getNodes())
    respT = np.zeros((samples, num_stations, C))
    arvR = np.zeros((samples, num_stations, C))
    util = np.zeros((samples, num_stations))
    for s in range(samples):
        solver = SolverSSA(model, samples=100000, seed=s + 1)
        respT[s, :, :] = solver.getAvgRespT()
        arvR[s, :, :] = solver.getAvgArvR()
        util[s, :] = np.sum(solver.getAvgUtil(), axis=1)
    return respT, arvR, util
