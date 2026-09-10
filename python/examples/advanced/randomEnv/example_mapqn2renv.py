"""
Example: Transform MMPP service queue to random environment model

This example shows how to transform a closed queueing network with MMPP2
service into a random environment model with exponential services modulated by
the MMPP phases.

MMPP2 (2-phase Markov Modulated Poisson Process) is a point process where
service completions occur at rates modulated by a 2-state Markov chain. This
transformation converts such a network into an equivalent random environment
model where:
  - Environment has 2 stages (one per MMPP phase)
  - Service rates are exponential with rates from MMPP D1 diagonal
  - Environment transitions follow MMPP D0 matrix structure
"""

from line_solver import (ClosedClass, CTMC, Delay, ENV, Exp, FLD,
                         GlobalConstants, JMT, MMPP2, MVA, Network, Queue,
                         SchedStrategy, VerboseLevel)
from line_solver.api.io.converters import mapqn2renv


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.SILENT)

    # Create closed queueing network with MMPP2 service
    model = Network('MMPP_ClosedQN')

    delay = Delay(model, 'Delay')                        # Think time station
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)    # Service station with MMPP2

    N = 5                                                # Number of jobs
    jobclass = ClosedClass(model, 'Class1', N, delay)
    delay.setService(jobclass, Exp(1.0))

    # Parameters: lambda0, lambda1, sigma01, sigma10
    #   Phase 0: slow service (rate 1.0)  -> high queue length
    #   Phase 1: fast service (rate 10.0) -> low queue length
    #   D0 = [[-1.2, 0.2], [0.3, -10.3]],  D1 = [[1.0, 0], [0, 10.0]]
    queue.setService(jobclass, MMPP2(1.0, 10.0, 0.2, 0.3))

    model.link(Network.serialRouting(delay, queue))

    print('=== MMPP2 Closed QN to Random Environment Transformation ===\n')
    print('Original Network:')
    print('  Topology: Delay -> MMPP2 Queue -> Delay (cyclic)')
    print('  Population: N = %d jobs' % N)
    print('  Think time: Exp(1.0)')
    print('  Service: MMPP2 with phases 0,1\n')

    print('Transforming to random environment model...')
    envModel = mapqn2renv(model)
    print('Transformation complete!\n')

    print('Environment Model Structure:')
    print('  Environment: MAPQN_Env with 2 stages')
    print('  Stage 0 (Phase 0): Closed QN with Exp service rate 1.0 (SLOW)')
    print('  Stage 1 (Phase 1): Closed QN with Exp service rate 10.0 (FAST)')
    print('  Transitions:')
    print('    Phase 0 -> Phase 1: Rate 0.2')
    print('    Phase 1 -> Phase 0: Rate 0.3\n')

    # Python keeps the stage list on the Environment itself; MATLAB reads the
    # same three columns off envGraph.Nodes.
    print('Environment stages:')
    for name, stage_type in zip(envModel._stage_names, envModel._stage_types):
        print('  %s: %s' % (name, stage_type))
    print()

    print('Stage network details:')
    for name, stage_model in zip(envModel._stage_names, envModel.getEnsemble()):
        if stage_model is not None:
            print('  %s: Network with %d nodes' % (name, stage_model.getNumberOfNodes()))
    print()

    # Solve the original MMPP2 model using JMT
    try:
        print('Solving original MMPP2 model using JMT...')
        print('Original MMPP2 Queue results:')
        print(JMT(model).getAvgTable())
    except Exception as err:
        print('JMT solver error: %s' % err)

    # Solve the environment model using ENV with FLD
    print('\nSolving environment model using ENV (with FLD)...')
    try:
        options = {'timespan': [0, float('inf')], 'iter_max': 100,
                   'iter_tol': 0.01, 'method': 'default', 'verbose': False}
        envSolver = ENV(envModel, lambda m: FLD(m, timespan=[0, 100], verbose=False),
                        options)
        QN, UN, _, TN = envSolver.getAvg()[:4]
        print('\nEnvironment Model AvgTable:')
        print(envSolver.getAvgTable())
        print('Note: The random environment model captures MMPP dynamics through')
        print('switching between phases with exponential service rates.\n')
    except Exception as err:
        print('ENV error: %s\n' % err)
        print('Note: ENV requires transient analysis support.')
        print('The environment model was created successfully.\n')

    # Solve the environment model using ENV with CTMC
    print('\nSolving environment model using ENV (with CTMC, cutoff=100)...')
    try:
        options = {'timespan': [0, float('inf')], 'iter_max': 100,
                   'iter_tol': 0.01, 'method': 'default', 'verbose': False}
        envSolverCTMC = ENV(envModel,
                            lambda m: CTMC(m, 'exact', timespan=[0, 100],
                                           cutoff=100, verbose=False),
                            options)
        QN_ctmc, UN_ctmc, _, TN_ctmc = envSolverCTMC.getAvg()[:4]
        print('\nEnvironment Model AvgTable (CTMC):')
        print(envSolverCTMC.getAvgTable())
    except Exception as err:
        print('ENV/CTMC error: %s\n' % err)

    # Solve individual stage networks using MVA (steady-state)
    print('Solving stage networks individually (steady-state with MVA)...')
    try:
        for name, stage_model in zip(envModel._stage_names, envModel.getEnsemble()):
            if stage_model is not None:
                print('\n  Stage %s:' % name)
                print(MVA(stage_model).getAvgTable())
    except Exception as err:
        print('Stage solver error: %s' % err)

    print('\n=== Transformation Complete ===')
