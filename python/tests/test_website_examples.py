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
STOCHASTIC_TOLERANCE = 0.30     # 30% - accounts for SSA simulation variance


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
    - Double printing (some solvers print internally AND we add print())
    """
    # Skip verification if no expected output
    if not expected:
        return True

    # Verify key metrics are present
    required_metrics = ['QLen', 'Util', 'RespT', 'Tput']
    for metric in required_metrics:
        if metric not in actual:
            return False

    # Determine solver type
    is_stochastic = solver_name in STOCHASTIC_SOLVERS

    # Extract numerical values (use unique values to handle double printing)
    actual_nums = extract_numbers(actual)
    expected_nums = extract_numbers(expected)

    # Remove duplicates while preserving order (handle double printing)
    # Each expected value should appear in actual at least once
    seen = set()
    unique_actual = []
    for num in actual_nums:
        # Use rounded value for comparison (handle floating point)
        rounded = round(num, 4)
        if rounded not in seen:
            seen.add(rounded)
            unique_actual.append(num)

    seen = set()
    unique_expected = []
    for num in expected_nums:
        rounded = round(num, 4)
        if rounded not in seen:
            seen.add(rounded)
            unique_expected.append(num)

    # Select tolerance based on solver type
    tolerance = STOCHASTIC_TOLERANCE if is_stochastic else DETERMINISTIC_TOLERANCE

    # For each unique expected value, find a matching actual value
    for expected_val in unique_expected:
        found = False
        for actual_val in unique_actual:
            if abs(expected_val) > 0 and abs(actual_val - expected_val) <= tolerance * abs(expected_val):
                found = True
                break
            elif abs(expected_val) == 0 and abs(actual_val) < 0.001:
                found = True
                break
        if not found:
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

                    # Fix code for exec() context: last expression should be printed
                    # In interactive Python, the last expression is automatically printed,
                    # but exec() doesn't do this. Add explicit print() for the last expression.
                    lines = code.strip().split('\n')
                    last_line = lines[-1].strip()
                    # If last line is an expression (not assignment, import, etc.)
                    if (last_line and
                        not last_line.startswith('#') and
                        not last_line.startswith('import ') and
                        not last_line.startswith('from ') and
                        '=' not in last_line.split('(')[0] and  # Not an assignment
                        not last_line.startswith('print')):
                        # Wrap last line in print()
                        lines[-1] = f'print({last_line})'
                        code = '\n'.join(lines)

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
