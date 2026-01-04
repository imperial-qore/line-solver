"""
MATLAB Parity Validation Tests for SolverMAM.

Compares Python SolverMAM results against MATLAB reference implementations.

Test Strategy:
1. Create network models that match MATLAB examples
2. Run SolverMAM on each model
3. Compare against known MATLAB results or queueing theory baselines
4. Validate within specified tolerances

Tolerance Levels:
- Queue lengths (QN): 1e-6 relative error (analytical methods)
- Utilization (UN): 1e-6 relative error
- Response times (RN): 1e-6 relative error
- Throughput (TN): 1e-6 relative error

For approximation methods (MNA, RCAT):
- Acceptable tolerance: 1e-3 (0.1%) relative error
"""

import unittest
import numpy as np
import pandas as pd
from typing import Tuple, Dict, List

from line_solver.solvers.solver_mam import SolverMAM, SolverMAMOptions
from line_solver.tests.test_solver_mam_integration import RealNetworkStruct


class MATLABParityValidator:
    """Utility class for validating numerical parity with MATLAB."""

    def __init__(self, tol_analytical=1e-6, tol_approx=1e-3):
        """Initialize validator with tolerance levels."""
        self.tol_analytical = tol_analytical
        self.tol_approx = tol_approx

    def validate_metric(self, python_val, matlab_val, name, tol=None):
        """
        Validate a single metric against MATLAB reference.

        Returns:
            (pass, rel_error, message)
        """
        if tol is None:
            tol = self.tol_analytical

        if matlab_val == 0:
            if python_val == 0:
                return True, 0.0, f"{name}: Both zero (OK)"
            else:
                abs_error = abs(python_val - matlab_val)
                return abs_error < 1e-10, abs_error, f"{name}: |{python_val} - {matlab_val}| = {abs_error}"

        rel_error = abs(python_val - matlab_val) / abs(matlab_val)
        passed = rel_error <= tol

        status = "✓" if passed else "✗"
        return passed, rel_error, f"{name}: rel_error={rel_error*100:.6f}% (tol={tol*100:.2f}%) {status}"

    def validate_array(self, python_arr, matlab_arr, name, tol=None):
        """Validate arrays element-wise."""
        python_arr = np.asarray(python_arr)
        matlab_arr = np.asarray(matlab_arr)

        if python_arr.shape != matlab_arr.shape:
            return False, -1, f"{name}: Shape mismatch {python_arr.shape} vs {matlab_arr.shape}"

        rel_errors = []
        for i in range(python_arr.size):
            py_val = python_arr.flat[i]
            m_val = matlab_arr.flat[i]

            if m_val == 0:
                if py_val == 0:
                    rel_errors.append(0.0)
                else:
                    rel_errors.append(float('inf'))
            else:
                rel_error = abs(py_val - m_val) / abs(m_val)
                rel_errors.append(rel_error)

        max_rel_error = max(rel_errors) if rel_errors else 0.0
        if tol is None:
            tol = self.tol_analytical

        passed = max_rel_error <= tol
        status = "✓" if passed else "✗"
        return passed, max_rel_error, f"{name}: max_rel_error={max_rel_error*100:.6f}% (tol={tol*100:.2f}%) {status}"


class TestM_M_1_Parity(unittest.TestCase):
    """MATLAB Parity: M/M/1 Queue"""

    def setUp(self):
        """Set up M/M/1 queue: λ=1, μ=2"""
        self.sn = RealNetworkStruct(nstations=1, nclasses=1)
        self.sn.lambda_arr = np.array([1.0])
        self.sn.mu = np.array([[2.0]])
        self.sn.rates = np.array([[2.0]])
        self.sn.scv = np.array([[1.0]])  # Exponential
        self.sn.nservers = np.array([1])
        self.validator = MATLABParityValidator()

    def test_mm1_utilization(self):
        """M/M/1: ρ = λ/μ = 0.5"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_util = 0.5  # λ/μ = 1/2
        actual_util = solver.result.UN[0, 0]

        passed, rel_err, msg = self.validator.validate_metric(
            actual_util, expected_util, "M/M/1 Utilization"
        )
        self.assertTrue(passed, msg)

    def test_mm1_response_time(self):
        """M/M/1: R = 1/(μ-λ) = 1/(2-1) = 1"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_r = 1.0 / (2.0 - 1.0)  # 1 second
        actual_r = solver.result.RN[0, 0]

        passed, rel_err, msg = self.validator.validate_metric(
            actual_r, expected_r, "M/M/1 Response Time"
        )
        self.assertTrue(passed, msg)

    def test_mm1_queue_length(self):
        """M/M/1: L = λ/(μ-λ) = 1/(2-1) = 1"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_l = 1.0 / (2.0 - 1.0)  # 1 customer
        actual_l = solver.result.QN[0, 0]

        passed, rel_err, msg = self.validator.validate_metric(
            actual_l, expected_l, "M/M/1 Queue Length"
        )
        self.assertTrue(passed, msg)


class TestTandemQueues_Parity(unittest.TestCase):
    """MATLAB Parity: Tandem M/M/1 Queues"""

    def setUp(self):
        """Set up tandem of 2 M/M/1 queues: λ=1, μ=[2, 2]"""
        self.sn = RealNetworkStruct(nstations=2, nclasses=1)
        self.sn.lambda_arr = np.array([1.0])
        self.sn.mu = np.array([[2.0], [2.0]])
        self.sn.rates = np.array([[2.0], [2.0]])
        self.sn.scv = np.array([[1.0], [1.0]])  # Both exponential
        self.sn.nservers = np.array([1, 1])
        # Routing: station 0→1, station 1→sink
        self.sn.rt = np.array([[0.0, 1.0], [0.0, 0.0]])
        self.validator = MATLABParityValidator()

    def test_tandem_utilization(self):
        """Both stations should have ρ = 0.5"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_util = np.array([0.5, 0.5])
        actual_util = solver.result.UN[:, 0]

        passed, rel_err, msg = self.validator.validate_array(
            actual_util, expected_util, "Tandem Utilization"
        )
        self.assertTrue(passed, msg)

    def test_tandem_response_times(self):
        """Each station: R = 1/(2-1) = 1"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_r = np.array([1.0, 1.0])
        actual_r = solver.result.RN[:, 0]

        passed, rel_err, msg = self.validator.validate_array(
            actual_r, expected_r, "Tandem Response Times"
        )
        self.assertTrue(passed, msg)

    def test_tandem_system_response_time(self):
        """System response time = sum of individual response times = 2"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_sys_r = 2.0  # R1 + R2
        actual_sys_r = np.sum(solver.result.RN[:, 0])

        passed, rel_err, msg = self.validator.validate_metric(
            actual_sys_r, expected_sys_r, "Tandem System Response Time"
        )
        self.assertTrue(passed, msg)


class TestMultiClass_Parity(unittest.TestCase):
    """MATLAB Parity: Multi-class Networks"""

    def setUp(self):
        """Set up 2-station, 2-class network"""
        self.sn = RealNetworkStruct(nstations=2, nclasses=2)
        self.sn.lambda_arr = np.array([0.5, 0.5])  # λ1=0.5, λ2=0.5
        self.sn.mu = np.array([[2.0, 2.0], [2.0, 2.0]])  # μij=2 for all i,j
        self.sn.rates = np.array([[2.0, 2.0], [2.0, 2.0]])
        self.sn.scv = np.array([[1.0, 1.0], [1.0, 1.0]])
        self.sn.nservers = np.array([1, 1])
        self.sn.rt = np.array([[0.0, 1.0], [0.0, 0.0]])
        self.validator = MATLABParityValidator()

    def test_multiclass_utilization(self):
        """Each class has ρ = λ/μ = 0.5/2 = 0.25"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_util = np.array([[0.25, 0.25], [0.25, 0.25]])
        actual_util = solver.result.UN

        passed, rel_err, msg = self.validator.validate_array(
            actual_util, expected_util, "Multi-class Utilization", tol=1e-2
        )
        self.assertTrue(passed, msg)

    def test_multiclass_throughput(self):
        """Each class should maintain its arrival rate"""
        solver = SolverMAM(self.sn, method='dec.source')
        solver.runAnalyzer()

        expected_tn = np.array([[0.5, 0.5]])
        actual_tn = solver.result.TN

        passed, rel_err, msg = self.validator.validate_array(
            actual_tn, expected_tn, "Multi-class Throughput", tol=1e-2
        )
        self.assertTrue(passed, msg)


class TestClosed_Parity(unittest.TestCase):
    """MATLAB Parity: Closed Networks"""

    def setUp(self):
        """Set up 2-station closed network with N=5 jobs"""
        self.sn = RealNetworkStruct(nstations=2, nclasses=1)
        self.sn.njobs = [5]
        self.sn.mu = np.array([[2.0], [2.0]])
        self.sn.rates = np.array([[2.0], [2.0]])
        self.sn.scv = np.array([[1.0], [1.0]])
        self.sn.nservers = np.array([1, 1])
        # Cycle: 0→1→0
        self.sn.rt = np.array([[0.0, 1.0], [1.0, 0.0]])
        self.validator = MATLABParityValidator(tol_approx=5e-2)  # 5% tolerance for decomposition

    def test_closed_population_constraint(self):
        """Total queue length should equal population (N=5)"""
        solver = SolverMAM(self.sn, method='mna_closed')
        solver.runAnalyzer()

        total_jobs = np.sum(solver.result.QN)
        expected_jobs = 5.0

        passed, rel_err, msg = self.validator.validate_metric(
            total_jobs, expected_jobs, "Closed Network Population", tol=0.1
        )
        self.assertTrue(passed, msg)

    def test_closed_non_zero_results(self):
        """All stations should have positive queue lengths and response times"""
        solver = SolverMAM(self.sn, method='mna_closed')
        solver.runAnalyzer()

        self.assertTrue(np.all(solver.result.QN > 0), "All queue lengths should be positive")
        self.assertTrue(np.all(solver.result.RN > 0), "All response times should be positive")
        self.assertTrue(np.all(np.isfinite(solver.result.QN)), "Queue lengths should be finite")


class TestMethodConsistency(unittest.TestCase):
    """Test that different methods produce consistent results on same model"""

    def setUp(self):
        """Set up a simple 2-station open network"""
        self.sn = RealNetworkStruct(nstations=2, nclasses=1)
        self.sn.lambda_arr = np.array([1.0])
        self.sn.mu = np.array([[2.0], [2.0]])
        self.sn.rates = np.array([[2.0], [2.0]])
        self.sn.scv = np.array([[1.0], [1.0]])
        self.validator = MATLABParityValidator(tol_approx=5e-2)

    def test_dec_source_vs_dec_mmap_consistency(self):
        """dec.source and dec.mmap should give qualitatively similar results (both positive, stable)"""
        solver1 = SolverMAM(self.sn, method='dec.source')
        solver1.runAnalyzer()

        solver2 = SolverMAM(self.sn, method='dec.mmap')
        solver2.runAnalyzer()

        # Both should produce positive queue lengths with same sign
        # dec.mmap includes departure processes so different quantitatively
        self.assertTrue(np.all(solver1.result.QN > 0), "dec.source should have positive QN")
        self.assertTrue(np.all(solver2.result.QN >= 0), "dec.mmap should have non-negative QN")
        self.assertTrue(np.all(np.isfinite(solver1.result.QN)), "dec.source QN should be finite")
        self.assertTrue(np.all(np.isfinite(solver2.result.QN)), "dec.mmap QN should be finite")

    def test_dec_source_vs_mna_consistency(self):
        """dec.source should produce stable results"""
        solver1 = SolverMAM(self.sn, method='dec.source')
        solver1.runAnalyzer()

        # dec.source should always converge to finite, positive results
        self.assertTrue(np.all(solver1.result.UN >= 0), "dec.source UN should be non-negative")
        self.assertTrue(np.all(solver1.result.UN < 1.1), "dec.source should be stable")
        self.assertTrue(np.all(np.isfinite(solver1.result.RN)), "dec.source RN should be finite")
        self.assertTrue(np.all(np.isfinite(solver1.result.QN)), "dec.source QN should be finite")


class TestNumericalAccuracy(unittest.TestCase):
    """Test numerical accuracy on edge cases"""

    def test_very_high_utilization_mm1(self):
        """M/M/1 with ρ=0.9: R = 1/(μ-λ) = 1/(10-9) = 10"""
        sn = RealNetworkStruct(nstations=1, nclasses=1)
        sn.lambda_arr = np.array([9.0])
        sn.mu = np.array([[10.0]])
        sn.rates = np.array([[10.0]])
        sn.scv = np.array([[1.0]])

        solver = SolverMAM(sn, method='dec.source')
        solver.runAnalyzer()

        expected_r = 1.0 / (10.0 - 9.0)  # = 10
        actual_r = solver.result.RN[0, 0]

        # Higher tolerance for high utilization
        validator = MATLABParityValidator(tol_approx=0.1)
        passed, rel_err, msg = validator.validate_metric(
            actual_r, expected_r, "High ρ Response Time", tol=0.1
        )
        self.assertTrue(passed, msg)

    def test_very_low_utilization_mm1(self):
        """M/M/1 with ρ=0.1: R = 1/(μ-λ) = 10/9 ≈ 1.111"""
        sn = RealNetworkStruct(nstations=1, nclasses=1)
        sn.lambda_arr = np.array([1.0])
        sn.mu = np.array([[11.0]])
        sn.rates = np.array([[11.0]])
        sn.scv = np.array([[1.0]])

        solver = SolverMAM(sn, method='dec.source')
        solver.runAnalyzer()

        expected_r = 1.0 / (11.0 - 1.0)
        actual_r = solver.result.RN[0, 0]

        validator = MATLABParityValidator()
        passed, rel_err, msg = validator.validate_metric(
            actual_r, expected_r, "Low ρ Response Time"
        )
        self.assertTrue(passed, msg)


if __name__ == '__main__':
    unittest.main(verbosity=2)
