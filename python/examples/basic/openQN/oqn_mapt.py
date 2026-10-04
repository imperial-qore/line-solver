"""
Open queueing network with a MAP_t arrival process and the Ko-Pender limits.

MAPt is a time-inhomogeneous Markovian arrival process: segment k covers
[breakpoints[k], breakpoints[k+1]) and carries the pair (D0[k], D1[k]), so the
stream is both non-renewal, through the modulating phase, and non-stationary,
through the schedule. Setting one phase recovers an NHPP; setting one segment
recovers an ordinary MAP.

SolverFLD's 'kp' method integrates the fluid and diffusion limits of Ko and
Pender, "Diffusion limits for the (MAP_t/Ph_t/inf)^N queueing network", Oper.
Res. Lett. 45 (2017) 248-253: the mean and the covariance of the queue length
are integrated jointly, so it is the only fluid method that returns a second
moment. For infinite-server stations the rate functions are affine in the state,
so both are exact rather than asymptotic -- Var(t) here is the exact variance of
the queue length, not an approximation.
"""

import numpy as np

from line_solver import (Delay, Exp, MAPt, Network, OpenClass, Sink, SolverFLD,
                         Source)

def oqn_mapt():
    model = Network('model')

    source = Source(model, 'Source')
    delay = Delay(model, 'Delay')
    sink = Sink(model, 'Sink')

    jobclass = OpenClass(model, 'OpenClass', 0)

    # Two segments of a 2-phase MAP, held for 1 and 1.5 time units and repeating.
    # The second segment runs the same phase graph at roughly twice the rate.
    breakpoints = [0.0, 1.0, 2.5]
    D0 = [np.array([[-5.0, 1.0], [2.0, -4.0]]),
          np.array([[-12.0, 3.0], [5.0, -9.0]])]
    D1 = [np.array([[3.0, 1.0], [1.0, 1.0]]),
          np.array([[7.0, 2.0], [2.0, 2.0]])]
    arrival = MAPt(breakpoints, D0, D1, True)

    source.setArrival(jobclass, arrival)
    delay.setService(jobclass, Exp(2.0))

    model.link(Network.serialRouting(source, delay, sink))

    # The arrival process rides along: the schedule is what the transient block
    # probes, and find_model searches a returned tuple element-wise, so handing
    # both back keeps the model exportable.
    return model, arrival


if __name__ == "__main__":
    # Built by the function above rather than inline: a model that exists only
    # under __main__ exposes nothing on import, so the JAVA and C++ parity rows
    # cannot export it and SKIP every solver -- a row that reads as coverage
    # while asserting nothing (see parity-static/_example_model_vendor.py).
    model, arrival = oqn_mapt()

    # Steady state of a cyclic schedule is the average over one period.
    # tol below the default 1e-4: the JAR's DormandPrince and MATLAB's ode15s
    # differ by ~0.2% at the default, which is integrator tolerance rather than a
    # modelling difference and would show up as a cross-codebase disagreement.
    print(SolverFLD(model, method='kp', tol=1e-9).getAvgTable())

    # Transient mean AND variance over two periods.
    solver = SolverFLD(model, method='kp', tol=1e-9)
    solver.options.timespan = (0.0, 5.0)
    solver.runAnalyzer()
    t = solver.result.t
    mean = solver.result.QNt[(1, 0)]
    var = solver.result.QVart[(1, 0)]
    print('\n   t     lambda(t)    QLen mean     QLen var')
    for probe in (0.5, 1.0, 1.5, 2.5, 3.0, 4.0, 5.0):
        print('%6.2f  %10.4f  %11.6f  %11.6f'
              % (probe, arrival.getRateAt(probe),
                 float(np.interp(probe, t, mean)),
                 float(np.interp(probe, t, var))))
