import shutil
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
import os
import sys
import unittest
import json
import warnings
import re
import numpy as np
from pathlib import Path
import subprocess

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# Install a kernel for the current Python interpreter to ensure notebooks run with the same environment
KERNEL_NAME = 'line_solver_test_kernel'
try:
    subprocess.run(
        [sys.executable, '-m', 'ipykernel', 'install', '--user', '--name', KERNEL_NAME, '--display-name', 'LINE Test Kernel'],
        check=True, capture_output=True
    )
except subprocess.CalledProcessError:
    # Fall back to default python3 kernel if installation fails
    KERNEL_NAME = 'python3'

TOLERANCE = 1e-6
RELATIVE_TOLERANCE = 1e-5


class TestGettingStarted(unittest.TestCase):

    notebooks_dir = 'examples/gettingstarted'
    regression_dir = 'tests/regression'

    def setUp(self):
        """Set up test environment."""
        self.working_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
        self.baseline_warnings = []

    def run_notebook(self, filename):
        os.chdir(self.working_dir)
        filepath = os.path.join(self.notebooks_dir, filename)
        print(f'Running notebook: {filename} Working directory: {self.working_dir}', flush=True)

        with open(filepath) as f:
            nb = nbformat.read(f, as_version=4)
            ep = ExecutePreprocessor(timeout=600, kernel_name=KERNEL_NAME)
            try:
                ep.preprocess(nb, {'metadata': {'path': os.path.dirname(filepath)}})
                print(f'Finished running notebook: {filename}', flush=True)
            except Exception as e:
                print(f"Notebook {filename} failed with error: {e}", flush=True)
            finally:
                if os.path.exists("_trial_temp"):
                    shutil.rmtree("_trial_temp")

        for cell in nb.cells:
            if cell.cell_type == 'code':
                for output in cell.get('outputs', []):
                    if output.output_type == 'error':
                        error_name = output.get("ename", "")
                        error_value = output.get("evalue", "")
                        self.fail(
                            f"Notebook {filename} failed.\n"
                            f"Error Name: {error_name}\n"
                            f"Error Value: {error_value}\n"
                        )

        self._compare_with_regression(filename, nb)

        return nb

    def _compare_with_regression(self, filename, nb):
        """Compare notebook results with saved regression data."""
        notebook_dir = os.path.dirname(filename)
        notebook_name = os.path.basename(filename).replace('.ipynb', '_regression.json')
        regression_path = os.path.join(self.working_dir, self.regression_dir, self.notebooks_dir, notebook_dir, notebook_name)

        if not os.path.exists(regression_path):
            warning_msg = f"No regression data found for {filename}. Run gen_all_regression.py to generate regression data."
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)
            return

        try:
            with open(regression_path, 'r') as f:
                regression_results = json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            warning_msg = f"Failed to load regression data for {filename}: {e}"
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)
            return

        current_results = self._extract_numerical_outputs(nb)

        differences = self._compare_results(regression_results, current_results, filename)

        if differences:
            warning_msg = f"Numerical differences detected in {filename}:\n" + "\n".join(differences)
            warnings.warn(warning_msg, UserWarning)
            self.baseline_warnings.append(warning_msg)

    def _extract_numerical_outputs(self, nb):
        """Extract numerical values from notebook outputs."""
        results = {}

        for cell_idx, cell in enumerate(nb.cells):
            if cell.cell_type == 'code':
                cell_results = []

                for output in cell.get('outputs', []):
                    if output.output_type in ['execute_result', 'display_data']:
                        if 'data' in output and 'text/plain' in output['data']:
                            text_output = output['data']['text/plain']
                            numerical_values = self._extract_numbers_from_text(text_output)
                            if numerical_values:
                                cell_results.extend(numerical_values)

                    elif output.output_type == 'stream' and output.name == 'stdout':
                        text_output = output['text']
                        numerical_values = self._extract_numbers_from_text(text_output)
                        if numerical_values:
                            cell_results.extend(numerical_values)

                if cell_results:
                    results[f'cell_{cell_idx}'] = cell_results

        return results

    def _extract_numbers_from_text(self, text):
        """Extract numerical values from text output using regex, filtering out runtimes and timestamps."""
        if not isinstance(text, str):
            return []

        text = self._filter_runtime_and_timestamps(text)

        number_pattern = r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?'

        matches = re.findall(number_pattern, text)

        numbers = []
        for match in matches:
            try:
                num = float(match)
                if not (np.isnan(num) or np.isinf(num)):
                    numbers.append(num)
            except (ValueError, OverflowError):
                continue

        return numbers

    def _filter_runtime_and_timestamps(self, text):
        """Filter out runtime and timestamp information from text."""
        if not isinstance(text, str):
            return text

        patterns_to_remove = [
            r'completed in \d+\.?\d*s\.?',
            r'\d{4}-\d{2}-\d{2}[\sT]\d{2}:\d{2}:\d{2}',
            r'\d{2}:\d{2}:\d{2}',
            r'\d{2}/\d{2}/\d{4}',
            r'\d{4}/\d{2}/\d{2}',
            r'/tmp/[^/\s]*\d{10,}[^/\s]*',
            r'/workspace/[^/\s]*\d{10,}[^/\s]*',
            r'JMT Model: [^\n]*',
            r'Execution time: \d+\.?\d*s?',
            r'Runtime: \d+\.?\d*s?',
            r'Time: \d+\.?\d*s?',
            r'\b\d{10,}\b',
        ]

        filtered_text = text
        for pattern in patterns_to_remove:
            filtered_text = re.sub(pattern, '', filtered_text, flags=re.IGNORECASE)

        return filtered_text

    def _compare_results(self, regression, current, filename):
        """Compare regression and current results, return list of differences."""
        differences = []

        for cell_key in regression:
            if cell_key not in current:
                differences.append(f"  {cell_key}: Missing in current results")
                continue

            regression_values = regression[cell_key]
            current_values = current[cell_key]

            if len(regression_values) != len(current_values):
                differences.append(
                    f"  {cell_key}: Different number of values "
                    f"(regression: {len(regression_values)}, current: {len(current_values)})"
                )
                continue

            for i, (regression_val, current_val) in enumerate(zip(regression_values, current_values)):
                abs_diff = abs(current_val - regression_val)
                rel_diff = abs_diff / abs(regression_val) if regression_val != 0 else abs_diff

                if abs_diff > TOLERANCE and rel_diff > RELATIVE_TOLERANCE:
                    differences.append(
                        f"  {cell_key}[{i}]: {current_val} vs {regression_val} "
                        f"(abs_diff: {abs_diff:.2e}, rel_diff: {rel_diff:.2e})"
                    )

        for cell_key in current:
            if cell_key not in regression:
                differences.append(f"  {cell_key}: Extra in current results")

        return differences

    def tearDown(self):
        """Print any baseline warnings at the end of each test."""
        if self.baseline_warnings:
            print("\n" + "="*60)
            print("REGRESSION COMPARISON WARNINGS:")
            print("="*60)
            for warning in self.baseline_warnings:
                print(warning)
            print("="*60)

    def test_tut01_mm1_basics(self):
        self.run_notebook("tut01_mm1_basics.ipynb")

    def test_tut02_mg1_multiclass_solvers(self):
        self.run_notebook("tut02_mg1_multiclass_solvers.ipynb")

    def test_tut03_repairmen(self):
        self.run_notebook("tut03_repairmen.ipynb")

    def test_tut04_lb_routing(self):
        self.run_notebook("tut04_lb_routing.ipynb")

    def test_tut05_completes_flag(self):
        self.run_notebook("tut05_completes_flag.ipynb")

    def test_tut06_cache_lru_zipf(self):
        self.run_notebook("tut06_cache_lru_zipf.ipynb")

    def test_tut07_respt_cdf(self):
        self.run_notebook("tut07_respt_cdf.ipynb")

    def test_tut08_opt_load_balancing(self):
        self.run_notebook("tut08_opt_load_balancing.ipynb")

    def test_tut09_dep_process_analysis(self):
        self.run_notebook("tut09_dep_process_analysis.ipynb")


if __name__ == '__main__':
    print(f"Running tests with tolerance settings:")
    print(f"  Absolute tolerance: {TOLERANCE}")
    print(f"  Relative tolerance: {RELATIVE_TOLERANCE}")
    print(f"  Regression directory: tests/regression/gettingstarted/")
    print("="*60)
    unittest.main(verbosity=2)