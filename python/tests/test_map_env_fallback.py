"""MAP/MMPP random-environment fallback of NetworkSolver.

A solver that cannot consume a non-renewal process solves the model through its
environment image instead of rejecting it. Targets are the MATLAB reference
values of the same models (matlab @NetworkSolver/mapEnvApprox.m), which the JAR
test jline/solvers/env/MapEnvFallbackTest.java also asserts.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Exp, MMPP2, MAP,
                         SchedStrategy, MVA, NC, CTMC, FLD)
from line_solver.api.io.converters import map2renv
from line_solver.api.sn import sn_map_modulation


def mmpp_closed():
    """Delay(Exp 1) -> Queue(MMPP2 1,10,.2,.3), N=5."""
    model = Network('mmppClosed')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Q1', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'C1', 5, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, MMPP2(1.0, 10.0, 0.2, 0.3))
    model.link(Network.serialRouting(delay, queue))
    return model


def gen_map_closed():
    """Closed network with a general MAP service (non-diagonal D1), N=4."""
    model = Network('genMap')
    delay = Delay(model, 'Think')
    queue = Queue(model, 'Q1', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'C1', 4, delay)
    delay.setService(jobclass, Exp(1.0))
    queue.setService(jobclass, MAP(np.array([[-3.0, 0.5], [0.2, -2.0]]),
                                   np.array([[2.0, 0.5], [1.0, 0.8]])))
    model.link(Network.serialRouting(delay, queue))
    return model


def test_modulation_records():
    mods = sn_map_modulation(mmpp_closed().getStruct())
    assert len(mods) == 1
    assert mods[0]['arrival'] is False
    assert mods[0]['classes'] == [0]
    assert mods[0]['order'] == 2
    assert bool(mods[0]['is_mmpp'])


def test_mmpp_service_intercepted_by_mva():
    solver = MVA(mmpp_closed())
    assert solver.needsMapEnv(solver.options)
    Q, U, R, T, A, W = solver.getAvg()
    # the 'auto' timescale test picks the rate-averaged limit here
    assert 'env.avg' in solver._result['method']
    assert T[0, 0] == pytest.approx(3.4427, abs=1e-3)
    assert Q[1, 0] == pytest.approx(1.5573, abs=1e-3)


def test_mmpp_service_intercepted_by_nc():
    solver = NC(mmpp_closed())
    assert solver.needsMapEnv(solver.options)
    Q, U, R, T, A, W = solver.getAvg()
    assert T[0, 0] == pytest.approx(3.4427, abs=1e-3)


def test_ctmc_solves_natively():
    solver = CTMC(mmpp_closed())
    assert not solver.needsMapEnv(solver.options)
    Q, U, R, T, A, W = solver.getAvg()
    assert T[0, 0] == pytest.approx(3.0952, abs=1e-3)


@pytest.mark.parametrize('method,expected', [('dec', 2.3424), ('avg', 3.4427)])
def test_forced_limits(method, expected):
    solver = MVA(mmpp_closed())
    solver.options.config = {'map_env_method': method}
    Q, U, R, T, A, W = solver.getAvg()
    assert T[0, 0] == pytest.approx(expected, abs=1e-3)


def test_map_env_off_restores_rejection():
    solver = MVA(mmpp_closed())
    solver.options.config = {'map_env': 'off'}
    assert not solver.needsMapEnv(solver.options)
    with pytest.raises(Exception):
        solver.getAvg()


def test_general_map_image_is_intensity_matched():
    model = gen_map_closed()
    env, info = map2renv(model)
    assert info['nstages'] == 2
    assert not info['is_mmpp']
    Q, U, R, T, A, W = MVA(model).getAvg()
    assert T[0, 0] == pytest.approx(1.9317, abs=1e-3)


def test_two_processes_give_product_stage_space():
    model = Network('twoMaps')
    delay = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    jobclass = ClosedClass(model, 'C1', 3, delay)
    delay.setService(jobclass, Exp(1.0))
    q1.setService(jobclass, MMPP2(1.0, 4.0, 0.5, 0.5))
    q2.setService(jobclass, MMPP2(2.0, 6.0, 0.4, 0.6))
    model.link(Network.serialRouting(delay, q1, q2))

    env, info = map2renv(model)
    assert info['nstages'] == 4
    assert info['orders'] == [2, 2]
    Q, U, R, T, A, W = MVA(model).getAvg()
    assert T[0, 0] == pytest.approx(1.5041, abs=1e-3)


def test_transient_capability_gates_meanfield():
    assert not MVA(mmpp_closed()).supportsTransientAnalysis()
    assert not NC(mmpp_closed()).supportsTransientAnalysis()
    assert CTMC(mmpp_closed()).supportsTransientAnalysis()
    assert FLD(mmpp_closed()).supportsTransientAnalysis()

    solver = MVA(mmpp_closed())
    solver.options.config = {'map_env_method': 'meanfield'}
    with pytest.raises(Exception):
        solver.getAvg()


def test_citations_report_the_environment_image():
    solver = MVA(mmpp_closed())
    solver.getAvg()
    keys = {e['key'] for e in solver.citations()}
    assert 'neut79' in keys
    assert 'CasT11' in keys
    assert 'yin.zhan98' in keys
