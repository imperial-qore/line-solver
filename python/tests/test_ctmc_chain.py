"""
Native-Python regression test: chain mode of SolverCTMC, i.e. the solver built
from a user-supplied MarkovProcess (CTMC) or MarkovChain (DTMC) instead of a
Network. The solver then skips state-space generation and solves the given
generator directly.

The reference is the closed-form stationary vector of the two-state chains used
here. The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_chain.m and the JAR twin is
SolverCTMCChainTest.java.
"""

import os
import sys

_WORKTREE_PY = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if sys.path and sys.path[0] != _WORKTREE_PY:
    sys.path.insert(0, _WORKTREE_PY)

import numpy as np
import pytest

from line_solver import MarkovChain, MarkovProcess, SolverCTMC

import line_solver as _ls
assert os.path.abspath(_ls.__file__).startswith(_WORKTREE_PY), (
    f'wrong line_solver imported: {_ls.__file__} (expected under {_WORKTREE_PY})')

TOL = 1e-9

# Two-state CTMC with pi = [2/7, 5/7]
Q = np.array([[-0.5, 0.5], [0.2, -0.2]])
# Two-state DTMC with pi = [0.8, 0.2]
P = np.array([[0.9, 0.1], [0.4, 0.6]])


def test_ctmc_chain_stationary_matches_closed_form():
    solver = SolverCTMC(MarkovProcess(Q))
    pi = solver.getProbSys()
    assert solver.isChainSolver()
    assert not solver.isDiscreteChain()
    assert pi[0] == pytest.approx(2.0 / 7.0, abs=TOL)
    assert pi[1] == pytest.approx(5.0 / 7.0, abs=TOL)


def test_ctmc_chain_returns_the_given_generator():
    solver = SolverCTMC(MarkovProcess(Q))
    infgen, event_filt = solver.getGenerator()
    assert np.allclose(infgen, Q, atol=TOL)
    # No state space attached, so the state space is the state indices
    space, local = solver.getStateSpace()
    assert space.shape == (2, 1)
    assert np.allclose(space.flatten(), [1, 2])
    assert len(local) == 1
    assert np.allclose(solver.getStateSpaceAggr(), space)


def test_ctmc_chain_prob_of_single_state():
    ctmc = MarkovProcess(Q, True, np.array([[0], [1]]))
    solver = SolverCTMC(ctmc)
    assert solver.getProb(1) == pytest.approx(5.0 / 7.0, abs=TOL)
    with pytest.raises(RuntimeError):
        solver.getProb(7)


# Cross-codebase reference: pi(1) of the CTMC and pi after 5 steps of the DTMC,
# both from the uniform start. MATLAB, JAR and Python assert these same numbers.
PI_CTMC_T1 = np.array([0.392125422241, 0.607874577759])
PI_DTMC_K5 = np.array([0.790625, 0.209375])


def test_chain_transient_matches_the_cross_codebase_reference():
    ctmc = SolverCTMC(MarkovProcess(Q))
    assert np.allclose(ctmc.getTranProbSys(1.0), PI_CTMC_T1, atol=1e-5)
    dtmc = SolverCTMC(MarkovChain(P))
    assert np.allclose(dtmc.getTranProbSys(5), PI_DTMC_K5, atol=1e-9)
    # the container methods answer the same distribution as the solver
    assert np.allclose(MarkovProcess(Q).transientProb(None, 1.0), PI_CTMC_T1, atol=1e-5)
    assert np.allclose(MarkovChain(P).transientProb(None, 5)[-1], PI_DTMC_K5, atol=1e-9)


def test_ctmc_chain_transient_converges_to_stationary():
    solver = SolverCTMC(MarkovProcess(Q))
    pi_t = solver.getTranProbSys(100.0)
    assert pi_t[0] == pytest.approx(2.0 / 7.0, abs=1e-3)
    assert pi_t[1] == pytest.approx(5.0 / 7.0, abs=1e-3)


def test_dtmc_chain_stationary_matches_closed_form():
    dtmc = MarkovChain(P)
    solver = SolverCTMC(dtmc)
    pi = solver.getProbSys()
    assert solver.isDiscreteChain()
    assert pi[0] == pytest.approx(0.8, abs=TOL)
    assert pi[1] == pytest.approx(0.2, abs=TOL)
    # MarkovChain.solve is the twin of MarkovProcess.solve
    pi_direct = np.asarray(dtmc.solve()).flatten()
    assert pi_direct[0] == pytest.approx(0.8, abs=TOL)
    assert pi_direct[1] == pytest.approx(0.2, abs=TOL)


def test_dtmc_chain_generator_is_the_uniformized_one():
    solver = SolverCTMC(MarkovChain(P))
    infgen, _ = solver.getGenerator()
    assert np.allclose(infgen, P - np.eye(2), atol=TOL)
    assert np.allclose(solver.getTransMat(), P, atol=TOL)


def test_dtmc_chain_transient_advances_one_step_per_unit_time():
    solver = SolverCTMC(MarkovChain(P))
    pi_t = solver.getTranProbSys(5)
    pi_k = np.array([0.5, 0.5])
    for _ in range(5):
        pi_k = pi_k @ P
    assert np.allclose(pi_t, pi_k, atol=TOL)
    with pytest.raises(RuntimeError):
        solver.getTranProbSys(2.5)


def test_chain_sample_stays_in_the_state_space():
    solver = SolverCTMC(MarkovProcess(Q), seed=1)
    sample = solver.sampleSys(50)
    assert sample.state.shape[0] == len(sample.t)
    assert np.all(np.isin(sample.state.flatten(), [1.0, 2.0]))
    assert np.all(np.diff(sample.t) > 0)


def test_markovprocess_api_methods():
    ctmc = MarkovProcess(Q)
    pi = np.array([2.0 / 7.0, 5.0 / 7.0])
    # transient endpoint and time average both converge to stationary
    assert np.allclose(ctmc.transient(None, 200.0), pi, atol=1e-6)
    time_avg, exit_dist = ctmc.timeAverage(None, 200.0)
    assert np.allclose(time_avg, pi, atol=1e-2)
    assert np.allclose(exit_dist, pi, atol=1e-6)
    # the constructor repairs the diagonal, so the carried generator is valid
    assert ctmc.isFeasible()
    # the embedded jump chain drops the holding times, so its stationary
    # vector differs from the CTMC one
    embedded = ctmc.toEmbedded()
    assert np.allclose(embedded.getTransMat(), np.array([[0.0, 1.0], [1.0, 0.0]]), atol=TOL)
    assert np.allclose(np.asarray(embedded.solve()).flatten(), [0.5, 0.5], atol=1e-9)
    # sensitivity of the stationary vector conserves total probability
    dQ = np.array([[-1.0, 1.0], [0.0, 0.0]])
    assert abs(np.asarray(ctmc.sens(dQ)).sum()) < 1e-8
    # stochastic complement of a single state is a 1x1 generator, and the
    # full form carries the blocks it was built from
    assert np.asarray(ctmc.stochComp(np.array([0]))).shape == (1, 1)
    blocks = ctmc.stochCompFull(np.array([0]))
    assert set(['S', 'Q11', 'Q12', 'Q21', 'Q22', 'T']) <= set(blocks.keys())
    assert np.allclose(blocks['S'], ctmc.stochComp(np.array([0])), atol=TOL)
    assert np.allclose(blocks['Q11'], [[-0.5]], atol=TOL)
    assert np.allclose(blocks['Q12'], [[0.5]], atol=TOL)
    assert np.allclose(blocks['Q21'], [[0.2]], atol=TOL)
    assert np.allclose(blocks['Q22'], [[-0.2]], atol=TOL)
    # T = Q12*inv(-Q22)*Q21 is the return path, and the complement is Q11 + T
    assert np.allclose(blocks['T'], blocks['Q12'] @ np.linalg.inv(-blocks['Q22']) @ blocks['Q21'], atol=TOL)
    assert np.allclose(blocks['S'], blocks['Q11'] + blocks['T'], atol=TOL)


def test_markovchain_api_methods():
    dtmc = MarkovChain(P)
    pi_t = dtmc.transient(np.array([0.5, 0.5]), 5)
    pi_k = np.array([0.5, 0.5])
    for _ in range(5):
        pi_k = pi_k @ P
    assert pi_t.shape == (6, 2)
    assert np.allclose(pi_t[-1], pi_k, atol=TOL)
    # from state 0 the target state 1 is hit after 1/0.1 steps on average
    h = dtmc.hittingTime([1])
    assert h[0] == pytest.approx(10.0, abs=1e-9)
    assert h[1] == pytest.approx(0.0, abs=TOL)
    assert dtmc.isFeasible()


def test_markovprocess_api_extras():
    ctmc = MarkovProcess(Q)
    # the two transient engines agree
    assert np.allclose(ctmc.transient(None, 1.0, 'foxglynn'),
                       ctmc.transient(None, 1.0, 'unif'), atol=1e-9)
    # the relative solution is the stationary vector scaled by p(refstate)
    rel = np.asarray(ctmc.solveRelative(0)).flatten()
    pi = np.asarray(ctmc.solve()).flatten()
    assert rel[0] == pytest.approx(1.0, abs=TOL)
    assert np.allclose(rel / rel.sum(), pi, atol=1e-9)
    # aggregation of a nearly-decomposable chain recovers the exact vector
    Qb = np.array([[-1.001, 1.0, 0.001, 0.0],
                   [1.0, -1.001, 0.0, 0.001],
                   [0.001, 0.0, -1.001, 1.0],
                   [0.0, 0.001, 1.0, -1.001]])
    block = MarkovProcess(Qb)
    exact = np.asarray(block.solve()).flatten()
    for method in ['courtois', 'kms', 'takahashi']:
        p_agg, eps, eps_max = block.aggregate([[0, 1], [2, 3]], method)
        assert np.allclose(np.asarray(p_agg).flatten(), exact, atol=1e-3)
        assert eps <= eps_max
    p_multi, _, _ = block.aggregate([[0, 1], [2, 3]], 'multi', [[0, 1]])
    assert np.allclose(np.asarray(p_multi).flatten(), exact, atol=1e-3)


def test_markovchain_api_extras():
    dtmc = MarkovChain(P)
    blocks = dtmc.stochCompFull([0])
    assert np.allclose(blocks['P11'], [[0.9]], atol=TOL)
    assert np.allclose(blocks['P12'], [[0.1]], atol=TOL)
    assert np.allclose(blocks['P21'], [[0.4]], atol=TOL)
    assert np.allclose(blocks['P22'], [[0.6]], atol=TOL)
    # S = P11 + P12*inv(I-P22)*P21, a 1x1 stochastic matrix here
    assert np.allclose(blocks['S'], blocks['P11'] + blocks['P12'] @ np.linalg.inv(
        np.eye(1) - blocks['P22']) @ blocks['P21'], atol=TOL)
    # read as the randomized image of a CTMC, the chain relaxes to its
    # stationary vector as t grows
    assert np.allclose(dtmc.transientUnif(None, 500.0), [0.8, 0.2], atol=1e-6)


def test_average_metrics_are_refused_in_chain_mode():
    solver = SolverCTMC(MarkovProcess(Q))
    with pytest.raises(RuntimeError):
        solver.getAvgTable()
    with pytest.raises(RuntimeError):
        solver.getAvgQLen()
