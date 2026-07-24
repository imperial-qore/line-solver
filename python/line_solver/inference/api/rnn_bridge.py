"""MATLAB-Python bridge for RNN demand estimation.

Called from MATLAB's estimator_rnn.m via system(). Reads input data from
a MAT file, trains the RNN, and writes results to an output MAT file.

Usage:
    python -m line_solver.inference.api.rnn_bridge <input.mat> <output.mat>

Input MAT file fields:
    traces     - 4D array (traceCount x S x M x R+1)
    numServers - 1D array (M,)

Output MAT file fields:
    demandEst  - 1D array (M,) of estimated mean service demands

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
This code is released under the 3-Clause BSD License.
"""

import sys
import numpy as np
import scipy.io as sio


def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.mat> <output.mat>", file=sys.stderr)
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    # Load input data
    data = sio.loadmat(input_path)
    traces = data['traces'].astype(np.float32)
    num_servers = data['numServers'].flatten().astype(np.float32)

    # Run RNN estimation
    from line_solver.inference.api._estimators import _rnn_data
    demand_est = _rnn_data(traces, num_servers)

    # Save results
    sio.savemat(output_path, {'demandEst': demand_est})


if __name__ == '__main__':
    main()
