"""Test suite for line_solver CLI functionality."""

import os
import sys
import unittest
import tempfile
import json
from io import StringIO
from pathlib import Path

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from line_solver.cli import (
    create_parser, detect_input_format, format_readable, format_json, format_csv,
    validate_solver_compatibility, JSIM_COMPATIBLE_SOLVERS, LQN_COMPATIBLE_SOLVERS,
    VALID_ANALYSIS_TYPES, ANALYSIS_SOLVER_COMPAT, ANALYSIS_REQUIRES_NODE,
    ANALYSIS_REQUIRES_CLASS, VALID_REWARD_NAMES,
    validate_analysis_types, validate_analysis_solver_compat, validate_analysis_params
)
from line_solver import Network, OpenClass, Source, Queue, Sink, Exp, SchedStrategy


class TestCLIParser(unittest.TestCase):
    """Test argument parser creation and parsing."""

    def setUp(self):
        """Set up test fixtures."""
        self.parser = create_parser()

    def test_parser_creation(self):
        """Test that parser is created successfully."""
        self.assertIsNotNone(self.parser)

    def test_help_parsing(self):
        """Test help parsing doesn't crash."""
        # This should not raise an exception
        try:
            self.parser.print_help()
        except SystemExit:
            # print_help calls sys.exit(0), which is expected
            pass

    def test_default_arguments(self):
        """Test default argument values."""
        args = self.parser.parse_args([])
        self.assertEqual(args.solver, 'auto')
        self.assertEqual(args.analysis, 'all')
        self.assertEqual(args.output, 'readable')
        self.assertEqual(args.verbosity, 'normal')
        self.assertIsNone(args.file)
        self.assertIsNone(args.port)

    def test_solver_argument(self):
        """Test solver argument parsing."""
        args = self.parser.parse_args(['-s', 'ctmc'])
        self.assertEqual(args.solver, 'ctmc')

    def test_analysis_argument(self):
        """Test analysis type argument parsing."""
        for analysis_type in ['all', 'avg', 'sys']:
            args = self.parser.parse_args(['-a', analysis_type])
            self.assertEqual(args.analysis, analysis_type)

    def test_output_format_argument(self):
        """Test output format argument parsing."""
        for fmt in ['readable', 'json', 'csv', 'pickle', 'mat']:
            args = self.parser.parse_args(['-o', fmt])
            self.assertEqual(args.output, fmt)


class TestFormatDetection(unittest.TestCase):
    """Test input format auto-detection."""

    def test_detect_jsim_formats(self):
        """Test detection of JSIM formats."""
        self.assertEqual(detect_input_format('model.jsim'), 'jsim')
        self.assertEqual(detect_input_format('model.jsimg'), 'jsimg')
        self.assertEqual(detect_input_format('model.jsimw'), 'jsimw')
        self.assertEqual(detect_input_format('model.JSIM'), 'jsim')

    def test_detect_lqn_formats(self):
        """Test detection of LQN formats."""
        self.assertEqual(detect_input_format('model.lqnx'), 'lqnx')
        self.assertEqual(detect_input_format('model.xml'), 'xml')
        self.assertEqual(detect_input_format('model.LQNX'), 'lqnx')

    def test_detect_other_formats(self):
        """Test detection of other formats."""
        self.assertEqual(detect_input_format('model.mat'), 'mat')
        self.assertEqual(detect_input_format('model.pkl'), 'pkl')

    def test_unsupported_format(self):
        """Test that unsupported formats raise ValueError."""
        with self.assertRaises(ValueError):
            detect_input_format('model.xyz')

    def test_none_filename(self):
        """Test that None filename returns None."""
        self.assertIsNone(detect_input_format(None))


class TestSolverCompatibility(unittest.TestCase):
    """Test solver compatibility validation."""

    def test_jsim_compatible_solvers(self):
        """Test that JSIM-compatible solvers pass validation."""
        for solver in JSIM_COMPATIBLE_SOLVERS:
            try:
                validate_solver_compatibility('jsimg', solver)
            except ValueError:
                self.fail(f"Solver {solver} should be compatible with JSIM format")

    def test_jsim_incompatible_solvers(self):
        """Test that non-JSIM solvers fail validation."""
        incompatible_solvers = ['ln', 'lqns', 'env']
        for solver in incompatible_solvers:
            with self.assertRaises(ValueError):
                validate_solver_compatibility('jsimg', solver)

    def test_lqn_compatible_solvers(self):
        """Test that LQN-compatible solvers pass validation."""
        for solver in LQN_COMPATIBLE_SOLVERS:
            try:
                validate_solver_compatibility('lqnx', solver)
            except ValueError:
                self.fail(f"Solver {solver} should be compatible with LQN format")

    def test_lqn_incompatible_solvers(self):
        """Test that non-LQN solvers fail validation."""
        incompatible_solvers = ['fluid', 'jmt', 'ssa']
        for solver in incompatible_solvers:
            with self.assertRaises(ValueError):
                validate_solver_compatibility('lqnx', solver)


class TestOutputFormatting(unittest.TestCase):
    """Test output formatting functions."""

    def setUp(self):
        """Create test data."""
        # Create a simple test model
        model = Network('TestModel')
        source = Source(model, 'Source')
        job_class = OpenClass(model, 'Class1')
        queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
        sink = Sink(model, 'Sink')

        source.setArrival(job_class, Exp(1.0))
        queue.setService(job_class, Exp(0.5))

        model.addLink(source, queue)
        model.addLink(queue, sink)

        # Run a quick analysis
        from line_solver import SolverMVA
        solver = SolverMVA(model)
        self.results = {
            'AvgTable': solver.avg_table(),
            'AvgSysTable': None
        }

    def test_readable_formatting(self):
        """Test readable format output."""
        output = format_readable(self.results)
        self.assertIsNotNone(output)
        self.assertIsInstance(output, str)
        self.assertIn('AVGTABLE', output)

    def test_json_formatting(self):
        """Test JSON format output."""
        output = format_json(self.results)
        self.assertIsNotNone(output)
        self.assertIsInstance(output, str)
        # Verify it's valid JSON
        parsed = json.loads(output)
        self.assertIn('AvgTable', parsed)

    def test_csv_formatting(self):
        """Test CSV format output."""
        output = format_csv(self.results)
        self.assertIsNotNone(output)
        self.assertIsInstance(output, str)
        self.assertIn('AvgTable', output)


class TestCLIWithModels(unittest.TestCase):
    """Test CLI functionality with actual models."""

    def test_simple_mm1_model(self):
        """Test CLI with a simple M/M/1 model."""
        from line_solver import SolverMVA

        # Create simple M/M/1 model
        model = Network('MM1')
        source = Source(model, 'Source')
        job_class = OpenClass(model, 'Class1')
        queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
        sink = Sink(model, 'Sink')

        source.setArrival(job_class, Exp(1.0))
        queue.setService(job_class, Exp(2.0))

        model.addLink(source, queue)
        model.addLink(queue, sink)

        # Test with MVA solver
        solver = SolverMVA(model)
        results = solver.avg_table()

        self.assertIsNotNone(results)
        # Check that results have expected columns
        self.assertIn('RespT', results.columns)
        self.assertIn('Util', results.columns)

    def test_model_export_import(self):
        """Test loading models from pickle format."""
        from line_solver import LINE
        import pickle

        # Note: LINE Network objects can't be pickled due to wrapper state,
        # so we just test that the load function can handle pickle format
        # by creating a simple dict and saving it

        test_data = {'model_name': 'test', 'description': 'test model'}
        with tempfile.NamedTemporaryFile(suffix='.pkl', delete=False) as f:
            pickle_file = f.name

        try:
            with open(pickle_file, 'wb') as f:
                pickle.dump(test_data, f)

            # Load it back using LINE.load
            loaded_data = LINE.load(pickle_file)
            self.assertIsNotNone(loaded_data)
            self.assertEqual(loaded_data, test_data)
        finally:
            if os.path.exists(pickle_file):
                os.unlink(pickle_file)


class TestExtendedAnalysisTypes(unittest.TestCase):
    """Test extended analysis type constants and validation."""

    def test_valid_analysis_types_defined(self):
        """Test that all expected analysis types are defined."""
        expected_types = [
            'all', 'avg', 'sys', 'stage', 'chain', 'node', 'nodechain',
            'cdf-respt', 'cdf-passt', 'perct-respt',
            'tran-avg', 'tran-cdf-respt', 'tran-cdf-passt',
            'prob', 'prob-aggr', 'prob-marg', 'prob-sys', 'prob-sys-aggr',
            'sample', 'sample-aggr', 'sample-sys', 'sample-sys-aggr',
            'reward', 'reward-steady', 'reward-value'
        ]
        for analysis_type in expected_types:
            self.assertIn(analysis_type, VALID_ANALYSIS_TYPES,
                         f"Analysis type {analysis_type} should be in VALID_ANALYSIS_TYPES")

    def test_validate_analysis_types_single(self):
        """Test validation of single analysis types."""
        for analysis_type in ['all', 'avg', 'sys', 'stage', 'cdf-respt', 'sample']:
            try:
                validate_analysis_types(analysis_type)
            except ValueError:
                self.fail(f"Analysis type {analysis_type} should be valid")

    def test_validate_analysis_types_comma_separated(self):
        """Test validation of comma-separated analysis types."""
        try:
            validate_analysis_types('avg,sys')
            validate_analysis_types('avg,stage,chain')
            validate_analysis_types('sample,sample-aggr,sample-sys')
        except ValueError:
            self.fail("Comma-separated analysis types should be valid")

    def test_validate_analysis_types_invalid(self):
        """Test that invalid analysis types are rejected."""
        with self.assertRaises(ValueError):
            validate_analysis_types('invalid-type')
        with self.assertRaises(ValueError):
            validate_analysis_types('avg,invalid')

    def test_analysis_requires_node_defined(self):
        """Test that node-requiring analysis types are defined."""
        expected = ['prob', 'prob-aggr', 'prob-marg', 'sample', 'sample-aggr']
        for analysis_type in expected:
            self.assertIn(analysis_type, ANALYSIS_REQUIRES_NODE,
                         f"{analysis_type} should require node index")

    def test_analysis_requires_class_defined(self):
        """Test that class-requiring analysis types are defined."""
        self.assertIn('prob-marg', ANALYSIS_REQUIRES_CLASS)


class TestAnalysisSolverCompatibility(unittest.TestCase):
    """Test analysis type and solver compatibility validation."""

    def test_sample_requires_ssa(self):
        """Test that sample analysis requires SSA solver."""
        # Should pass for ssa
        try:
            validate_analysis_solver_compat(['sample'], 'ssa')
        except ValueError:
            self.fail("sample analysis should be compatible with ssa solver")

        # Should fail for non-ssa solvers
        with self.assertRaises(ValueError):
            validate_analysis_solver_compat(['sample'], 'mva')
        with self.assertRaises(ValueError):
            validate_analysis_solver_compat(['sample'], 'fluid')

    def test_reward_requires_ctmc(self):
        """Test that reward analysis requires CTMC solver."""
        # Should pass for ctmc
        try:
            validate_analysis_solver_compat(['reward'], 'ctmc')
            validate_analysis_solver_compat(['reward-steady'], 'ctmc')
        except ValueError:
            self.fail("reward analysis should be compatible with ctmc solver")

        # Should fail for non-ctmc solvers
        with self.assertRaises(ValueError):
            validate_analysis_solver_compat(['reward'], 'mva')
        with self.assertRaises(ValueError):
            validate_analysis_solver_compat(['reward-steady'], 'ssa')

    def test_perct_respt_requires_mam(self):
        """Test that perct-respt analysis requires MAM solver."""
        try:
            validate_analysis_solver_compat(['perct-respt'], 'mam')
        except ValueError:
            self.fail("perct-respt analysis should be compatible with mam solver")

        with self.assertRaises(ValueError):
            validate_analysis_solver_compat(['perct-respt'], 'mva')

    def test_avg_compatible_with_all_solvers(self):
        """Test that avg analysis is compatible with common solvers."""
        for solver in ['mva', 'ctmc', 'fluid', 'ssa', 'nc']:
            try:
                validate_analysis_solver_compat(['avg'], solver)
            except ValueError:
                self.fail(f"avg analysis should be compatible with {solver}")


class TestAnalysisParamValidation(unittest.TestCase):
    """Test analysis parameter validation."""

    def test_node_required_for_sample(self):
        """Test that sample analysis requires node index."""
        with self.assertRaises(ValueError):
            validate_analysis_params(['sample'], None, None, None)

    def test_node_provided_for_sample(self):
        """Test that sample analysis passes with node index."""
        try:
            validate_analysis_params(['sample'], 1, None, None)
        except ValueError:
            self.fail("sample analysis should pass with node index provided")

    def test_class_required_for_prob_marg(self):
        """Test that prob-marg analysis requires class index."""
        with self.assertRaises(ValueError):
            validate_analysis_params(['prob-marg'], 1, None, None)

    def test_class_provided_for_prob_marg(self):
        """Test that prob-marg analysis passes with class index."""
        try:
            validate_analysis_params(['prob-marg'], 1, 0, None)
        except ValueError:
            self.fail("prob-marg analysis should pass with node and class index provided")

    def test_reward_name_valid(self):
        """Test that valid reward names are accepted."""
        for name in VALID_REWARD_NAMES:
            try:
                validate_analysis_params(['reward-value'], None, None, name)
            except ValueError:
                self.fail(f"Reward name {name} should be valid")

    def test_avg_no_params_required(self):
        """Test that avg analysis doesn't require params."""
        try:
            validate_analysis_params(['avg'], None, None, None)
        except ValueError:
            self.fail("avg analysis should not require any params")


class TestNewParserArguments(unittest.TestCase):
    """Test new CLI argument parsing."""

    def setUp(self):
        """Set up test fixtures."""
        self.parser = create_parser()

    def test_node_argument(self):
        """Test node index argument parsing."""
        args = self.parser.parse_args(['--node', '1'])
        self.assertEqual(args.node, 1)

    def test_class_idx_argument(self):
        """Test class index argument parsing."""
        args = self.parser.parse_args(['--class-idx', '0'])
        self.assertEqual(args.class_idx, 0)

    def test_events_argument(self):
        """Test events argument parsing."""
        args = self.parser.parse_args(['--events', '5000'])
        self.assertEqual(args.events, 5000)

    def test_events_default(self):
        """Test events default value."""
        args = self.parser.parse_args([])
        self.assertEqual(args.events, 1000)

    def test_percentiles_argument(self):
        """Test percentiles argument parsing."""
        args = self.parser.parse_args(['--percentiles', '50,90,95,99'])
        self.assertEqual(args.percentiles, '50,90,95,99')

    def test_percentiles_default(self):
        """Test percentiles default value."""
        args = self.parser.parse_args([])
        self.assertEqual(args.percentiles, '50,90,95,99')

    def test_reward_name_argument(self):
        """Test reward name argument parsing."""
        args = self.parser.parse_args(['--reward-name', 'QLen'])
        self.assertEqual(args.reward_name, 'QLen')

    def test_extended_analysis_types(self):
        """Test parsing extended analysis types."""
        for analysis_type in ['stage', 'chain', 'cdf-respt', 'sample', 'reward']:
            args = self.parser.parse_args(['-a', analysis_type])
            self.assertEqual(args.analysis, analysis_type)

    def test_multi_analysis_parsing(self):
        """Test comma-separated analysis types."""
        args = self.parser.parse_args(['-a', 'avg,sys,stage'])
        self.assertEqual(args.analysis, 'avg,sys,stage')

    def test_mam_solver(self):
        """Test MAM solver argument."""
        args = self.parser.parse_args(['-s', 'mam'])
        self.assertEqual(args.solver, 'mam')


def run_tests():
    """Run all tests."""
    # Create test suite
    loader = unittest.TestLoader()
    suite = unittest.TestSuite()

    # Add tests
    suite.addTests(loader.loadTestsFromTestCase(TestCLIParser))
    suite.addTests(loader.loadTestsFromTestCase(TestFormatDetection))
    suite.addTests(loader.loadTestsFromTestCase(TestSolverCompatibility))
    suite.addTests(loader.loadTestsFromTestCase(TestOutputFormatting))
    suite.addTests(loader.loadTestsFromTestCase(TestCLIWithModels))
    suite.addTests(loader.loadTestsFromTestCase(TestExtendedAnalysisTypes))
    suite.addTests(loader.loadTestsFromTestCase(TestAnalysisSolverCompatibility))
    suite.addTests(loader.loadTestsFromTestCase(TestAnalysisParamValidation))
    suite.addTests(loader.loadTestsFromTestCase(TestNewParserArguments))

    # Run tests
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
