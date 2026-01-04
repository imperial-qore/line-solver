"""
Test that website examples produce correct outputs.

This test verifies that all examples shown on the LINE website solver*.html
pages produce the expected outputs. The test reads examples from the
website_examples.json file generated from the documentation HTML pages.
"""

import json
import re
import sys
import io
from pathlib import Path
import pytest


# Stochastic solvers require looser tolerance
STOCHASTIC_SOLVERS = {'JMT', 'SSA', 'DES'}
DETERMINISTIC_TOLERANCE = 0.01  # 1%
STOCHASTIC_TOLERANCE = 0.50     # 50% - high tolerance for simulation variance with finite samples


def find_json_file():
    """Find the website_examples.json file."""
    # Try from python/tests/ directory
    path1 = Path(__file__).parent.parent.parent / 'doc' / 'website_examples.json'
    if path1.exists():
        return path1

    # Try from python/ directory
    path2 = Path(__file__).parent.parent / 'doc' / 'website_examples.json'
    if path2.exists():
        return path2

    raise FileNotFoundError("Cannot find website_examples.json")


def extract_table_lines(text):
    """
    Extract only table data rows from output text.

    Table rows start with an index number followed by station/class data.
    This filters out progress lines (like "SSA samples: 100 200 ...") and
    informational lines (like "JMT Model: /tmp/...").
    """
    table_lines = []
    for line in text.split('\n'):
        stripped = line.strip()
        # Table data rows start with a digit followed by space and station name
        # e.g., "0  Source   Class1  0.0000  ..."
        if stripped and stripped[0].isdigit():
            # Check if this looks like a table row (has multiple columns)
            parts = stripped.split()
            if len(parts) >= 4 and not stripped.startswith('SSA') and not stripped.startswith('DES'):
                table_lines.append(stripped)
    return '\n'.join(table_lines)


def extract_numbers(text):
    """Extract all numerical values from text."""
    numbers = []
    pattern = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'
    matches = re.findall(pattern, text)

    for match in matches:
        try:
            num = float(match)
            if num != 0.0 and not (num != num):  # Skip zeros and NaN
                numbers.append(num)
        except ValueError:
            continue

    return numbers


def verify_output(actual, expected, solver_name):
    """
    Verify that actual output matches expected output.

    This function handles:
    - Stochastic solvers (JMT, SSA, DES) with tolerance for numerical values
    - Exact solvers (MVA, CTMC, etc.) with exact numerical matching
    - Whitespace and formatting variations
    """
    # Skip verification if no expected output
    if not expected:
        return True

    # Determine solver type
    is_stochastic = solver_name in STOCHASTIC_SOLVERS

    # Extract only table data rows to ignore progress/info lines
    actual_table = extract_table_lines(actual)
    expected_table = extract_table_lines(expected)

    # Extract numerical values from table rows only
    actual_nums = extract_numbers(actual_table)
    expected_nums = extract_numbers(expected_table)

    # Check if we have the same number of metrics
    if len(actual_nums) != len(expected_nums):
        return False

    # Select tolerance based on solver type
    tolerance = STOCHASTIC_TOLERANCE if is_stochastic else DETERMINISTIC_TOLERANCE

    # Compare numerical values
    for actual_val, expected_val in zip(actual_nums, expected_nums):
        if abs(actual_val - expected_val) > tolerance * abs(expected_val):
            return False

    # Verify key metrics are present
    required_metrics = ['QLen', 'Util', 'RespT', 'Tput']
    for metric in required_metrics:
        if metric not in actual:
            return False

    return True


def load_examples():
    """Load examples from JSON file."""
    json_path = find_json_file()
    with open(json_path, 'r', encoding='utf-8') as f:
        return json.load(f)


def extract_first_line(text):
    """Get the first non-empty line."""
    lines = text.split('\n')
    for line in lines:
        line = line.strip()
        if line:
            return line
    return ''


class TestWebsiteExamples:
    """Test class for website examples."""

    @classmethod
    def setup_class(cls):
        """Load examples once for all tests."""
        cls.examples = load_examples()

    def test_all_examples(self):
        """Test all Python examples from the website."""
        total_tests = 0
        passed_tests = 0
        failed_tests = []

        print("\n========================================")
        print("Testing Website Examples (Python)")
        print("========================================\n")

        for solver_name, examples in self.examples.items():
            # Find Python examples
            python_examples = [ex for ex in examples if ex['lang'] == 'python']

            if not python_examples:
                continue

            print(f"Testing {solver_name} solver ({len(python_examples)} examples)...")

            for i, example in enumerate(python_examples, 1):
                test_name = f"{solver_name}_python_{i}"
                total_tests += 1

                try:
                    # Execute the example code
                    code = example['code']
                    expected_output = example['output']

                    # Capture output
                    old_stdout = sys.stdout
                    sys.stdout = captured_output = io.StringIO()

                    try:
                        # Execute code
                        exec_globals = {}
                        exec(code, exec_globals)
                    finally:
                        sys.stdout = old_stdout

                    actual_output = captured_output.getvalue()

                    # Verify output
                    if verify_output(actual_output, expected_output, solver_name):
                        passed_tests += 1
                        print(f"  ✓ {test_name}: PASSED")
                    else:
                        failed_tests.append(test_name)
                        print(f"  ✗ {test_name}: FAILED (output mismatch)")
                        if expected_output:
                            print(f"    Expected output snippet: {extract_first_line(expected_output)}")

                except Exception as e:
                    failed_tests.append(test_name)
                    print(f"  ✗ {test_name}: FAILED ({str(e)})")

            print()

        # Print summary
        print("========================================")
        print("Summary:")
        print(f"  Total tests: {total_tests}")
        print(f"  Passed: {passed_tests}")
        print(f"  Failed: {total_tests - passed_tests}")

        if failed_tests:
            print("\nFailed tests:")
            for test in failed_tests:
                print(f"  - {test}")

        print("========================================\n")

        # Assert all tests passed
        assert passed_tests == total_tests, f"{total_tests - passed_tests} out of {total_tests} tests failed."


if __name__ == '__main__':
    # Run tests directly
    test = TestWebsiteExamples()
    test.setup_class()
    test.test_all_examples()
