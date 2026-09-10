"""
One-line Tandem Open Queueing Network

This example shows a compact solution of a tandem open queueing network.
"""
from line_solver import *
import numpy as np

if __name__ == "__main__":
    print('This example shows a compact solution of a tandem open queueing network.')

    # Model parameters
    lambda_rates = np.array([1, 2]) / 50  # lambda(r) - arrival rate of class r
    D = np.array([[10, 5], [5, 9]])  # D(i,r) - mean demand of class r at station i
    Z = np.array([91, 92])  # Z(r) - mean service time of class r at delay station

    # Create tandem PS network with infinite servers
    model = Network.tandem_ps_inf(lambda_rates, D, Z)

    # Solve with MVA
    solver = MVA(model)
    avg_table = solver.avg_table()
    print(avg_table)

    # 1-line version: avg_table = MVA(Network.tandem_ps_inf(lambda_rates, D, Z)).avg_table()
