#!/usr/bin/env python3
"""
Example 10: Basic layered queueing network
This example demonstrates a simple client-server application with two tiers
"""

import sys
import os

# Add the line_solver package to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from line_solver.layered import LayeredNetwork, Processor, Task, Entry, Activity
from line_solver.constants import SchedStrategy  
from line_solver.distributions import Exp
from line_solver.solvers import LN, MVA

def main():
    # Create the layered network model
    model = LayeredNetwork('ClientDBSystem')

    # Create processors
    P1 = Processor(model, 'ClientProcessor', 1, SchedStrategy.PS)
    P2 = Processor(model, 'DBProcessor', 1, SchedStrategy.PS)

    # Create tasks
    T1 = Task(model, 'ClientTask', 10, SchedStrategy.REF).on(P1)
    T1.set_think_time(Exp.fit_mean(5.0))  # 5-second think time
    T2 = Task(model, 'DBTask', float('inf'), SchedStrategy.INF).on(P2)

    # Create entries that represent service interfaces
    E1 = Entry(model, 'ClientEntry').on(T1)
    E2 = Entry(model, 'DBEntry').on(T2)

    # Define activities that specify the work performed and synchronous calls
    # Client activity: processes request and calls DB
    A1 = Activity(model, 'ClientActivity', Exp.fit_mean(1.0)).on(T1)
    A1.bound_to(E1).synch_call(E2, 2.5)  # 2.5 DB calls on average

    # DB activity: processes database request
    A2 = Activity(model, 'DBActivity', Exp.fit_mean(0.8)).on(T2)
    A2.bound_to(E2).replies_to(E2)

    # Solve the layered network using the LN solver with MVA applied to each layer
    solver = LN(model, lambda m: MVA(m))
    avg_table = solver.avg_table()
    print(avg_table)

if __name__ == '__main__':
    main()