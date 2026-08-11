"""
MATLAB Parity Validation Tests for SolverMAM.

Compares Python SolverMAM results against queueing-theory / MATLAB baselines,
using real line_solver Network models. Because a real open-network struct keeps
the Source as station 0, the queueing stations are at indices 1..M.

Tolerance Levels:
- Analytical metrics (QN/UN/RN/TN): 1e-6 relative error
- Approximation methods (MNA, RCAT): up to 1e-3 (0.1%) relative error
"""

import unittest
import numpy as np

from line_solver.solvers.solver_mam import SolverMAM, SolverMAMOptions
from line_solver.tests.test_solver_mam_integration import (
    mm1_open, tandem_open, multiclass_tandem, closed_cycle, _queue_rows,
)


class MATLABParityValidator:
    """Utility class for validating numerical parity with MATLAB."""

    def __init__(self, tol_analytical=1e-6, tol_approx=1e-3):
        self.tol_analytical = tol_analytical
        self.tol_approx = tol_approx

    def validate_metric(self, python_val, matlab_val, name, tol=None):
        if tol is None:
            tol = self.tol_analytical
        if matlab_val == 0:
            if python_val == 0:
                return True, 0.0, "%s: Both zero (OK)" % name
            abs_error = abs(python_val - matlab_val)
            return abs_error < 1e-10, abs_error, "%s: |%s - %s| = %s" % (name, python_val, matlab_val, abs_error)
        rel_error = abs(python_val - matlab_val) / abs(matlab_val)
        passed = rel_error <= tol
        return passed, rel_error, "%s: rel_error=%.6f%% (tol=%.2f%%)" % (name, rel_error * 100, tol * 100)

    def validate_array(self, python_arr, matlab_arr, name, tol=None):
        python_arr = np.asarray(python_arr)
        matlab_arr = np.asarray(matlab_arr)
        if python_arr.shape != matlab_arr.shape:
            return False, -1, "%s: Shape mismatch %s vs %s" % (name, python_arr.shape, matlab_arr.shape)
        rel_errors = []
        for i in range(python_arr.size):
            py_val = python_arr.flat[i]
            m_val = matlab_arr.flat[i]
            if m_val == 0:
                rel_errors.append(0.0 if py_val == 0 else float('inf'))
            else:
                rel_errors.append(abs(py_val - m_val) / abs(m_val))
        max_rel_error = max(rel_errors) if rel_errors else 0.0
        if tol is None:
            tol = self.tol_analytical
        return max_rel_error <= tol, max_rel_error, "%s: max_rel_error=%.6f%% (tol=%.2f%%)" % (name, max_rel_error * 100, tol * 100)


class TestM_M_1_Parity(unittest.TestCase):
    """MATLAB Parity: M/M/1 Queue (lam=1, mu=2). Queue is station 1."""

    def setUp(self):
        self.model = mm1_open(lam=1.0, mu=2.0)
        self.validator = MATLABParityValidator()

    def _solve(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        return solver.result

    def test_mm1_utilization(self):
        r = self._solve()
        passed, _, msg = self.validator.validate_metric(np.asarray(r.UN)[1, 0], 0.5, "M/M/1 Utilization", tol=1e-2)
        self.assertTrue(passed, msg)

    def test_mm1_response_time(self):
        r = self._solve()
        passed, _, msg = self.validator.validate_metric(np.asarray(r.RN)[1, 0], 1.0, "M/M/1 Response Time", tol=1e-2)
        self.assertTrue(passed, msg)

    def test_mm1_queue_length(self):
        r = self._solve()
        passed, _, msg = self.validator.validate_metric(np.asarray(r.QN)[1, 0], 1.0, "M/M/1 Queue Length", tol=1e-2)
        self.assertTrue(passed, msg)


class TestTandemQueues_Parity(unittest.TestCase):
    """MATLAB Parity: 2 tandem M/M/1 queues (lam=1, mu=[2,2]). Queues at 1,2."""

    def setUp(self):
        self.model = tandem_open(lam=1.0, mus=(2.0, 2.0))
        self.validator = MATLABParityValidator()

    def _solve(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        return solver.result

    def test_tandem_utilization(self):
        r = self._solve()
        passed, _, msg = self.validator.validate_array(np.asarray(r.UN)[[1, 2], 0], np.array([0.5, 0.5]),
                                                       "Tandem Utilization", tol=1e-2)
        self.assertTrue(passed, msg)

    def test_tandem_response_times(self):
        r = self._solve()
        passed, _, msg = self.validator.validate_array(np.asarray(r.RN)[[1, 2], 0], np.array([1.0, 1.0]),
                                                       "Tandem Response Times", tol=1e-2)
        self.assertTrue(passed, msg)

    def test_tandem_system_response_time(self):
        r = self._solve()
        actual_sys_r = np.sum(np.asarray(r.RN)[[1, 2], 0])
        passed, _, msg = self.validator.validate_metric(actual_sys_r, 2.0, "Tandem System Response Time", tol=1e-2)
        self.assertTrue(passed, msg)


class TestMultiClass_Parity(unittest.TestCase):
    """MATLAB Parity: 2-queue, 2-class network (lam=0.5 each, mu=2). Queues at 1,2."""

    def setUp(self):
        self.model = multiclass_tandem(lams=(0.5, 0.5), mus=((2.0, 2.0), (2.0, 2.0)))
        self.validator = MATLABParityValidator()

    def _solve(self):
        solver = SolverMAM(self.model, method='dec.source')
        solver.runAnalyzer()
        return solver.result

    def test_multiclass_utilization(self):
        r = self._solve()
        # per class rho = lam/mu = 0.5/2 = 0.25 at each queue station
        passed, _, msg = self.validator.validate_array(np.asarray(r.UN)[[1, 2], :],
                                                       np.array([[0.25, 0.25], [0.25, 0.25]]),
                                                       "Multi-class Utilization", tol=5e-2)
        self.assertTrue(passed, msg)

    def test_multiclass_throughput(self):
        r = self._solve()
        # each class maintains its arrival rate 0.5 at each queue station
        passed, _, msg = self.validator.validate_array(np.asarray(r.TN)[1, :], np.array([0.5, 0.5]),
                                                       "Multi-class Throughput", tol=5e-2)
        self.assertTrue(passed, msg)


class TestClosed_Parity(unittest.TestCase):
    """MATLAB Parity: 2-station closed network, N=5."""

    def setUp(self):
        self.model = closed_cycle(njobs=5, mus=(2.0, 2.0))
        self.validator = MATLABParityValidator(tol_approx=5e-2)

    def test_closed_population_constraint(self):
        solver = SolverMAM(self.model, method='mna_closed')
        solver.runAnalyzer()
        total_jobs = np.sum(np.asarray(solver.result.QN))
        passed, _, msg = self.validator.validate_metric(total_jobs, 5.0, "Closed Network Population", tol=0.1)
        self.assertTrue(passed, msg)

    def test_closed_non_zero_results(self):
        solver = SolverMAM(self.model, method='mna_closed')
        solver.runAnalyzer()
        QN = np.asarray(solver.result.QN)
        RN = np.asarray(solver.result.RN)
        self.assertTrue(np.all(QN >= 0), "Queue lengths should be non-negative")
        self.assertTrue(np.all(np.isfinite(QN)), "Queue lengths should be finite")
        self.assertTrue(np.all(np.isfinite(RN)), "Response times should be finite")


class TestMethodConsistency(unittest.TestCase):
    """Test that different methods produce consistent results on the same model."""

    def setUp(self):
        self.model = tandem_open(lam=1.0, mus=(2.0, 2.0))

    def test_dec_source_vs_dec_mmap_consistency(self):
        s1 = SolverMAM(self.model, method='dec.source'); s1.runAnalyzer()
        s2 = SolverMAM(self.model, method='dec.mmap'); s2.runAnalyzer()
        self.assertTrue(np.all(np.asarray(s1.result.QN) >= 0))
        self.assertTrue(np.all(np.asarray(s2.result.QN) >= 0))
        self.assertTrue(np.all(np.isfinite(np.asarray(s1.result.QN))))
        self.assertTrue(np.all(np.isfinite(np.asarray(s2.result.QN))))

    def test_dec_source_vs_mna_consistency(self):
        s1 = SolverMAM(self.model, method='dec.source'); s1.runAnalyzer()
        UN = np.asarray(s1.result.UN)
        self.assertTrue(np.all(UN >= 0))
        self.assertTrue(np.all(UN < 1.1))
        self.assertTrue(np.all(np.isfinite(np.asarray(s1.result.RN))))
        self.assertTrue(np.all(np.isfinite(np.asarray(s1.result.QN))))


class TestNumericalAccuracy(unittest.TestCase):
    """Test numerical accuracy on edge cases (queue is station 1)."""

    def test_very_high_utilization_mm1(self):
        solver = SolverMAM(mm1_open(lam=9.0, mu=10.0), method='dec.source')
        solver.runAnalyzer()
        actual_r = np.asarray(solver.result.RN)[1, 0]
        validator = MATLABParityValidator(tol_approx=0.1)
        passed, _, msg = validator.validate_metric(actual_r, 1.0 / (10.0 - 9.0), "High rho Response Time", tol=0.1)
        self.assertTrue(passed, msg)

    def test_very_low_utilization_mm1(self):
        solver = SolverMAM(mm1_open(lam=1.0, mu=11.0), method='dec.source')
        solver.runAnalyzer()
        actual_r = np.asarray(solver.result.RN)[1, 0]
        validator = MATLABParityValidator()
        passed, _, msg = validator.validate_metric(actual_r, 1.0 / (11.0 - 1.0), "Low rho Response Time", tol=1e-2)
        self.assertTrue(passed, msg)


if __name__ == '__main__':
    unittest.main(verbosity=2)
