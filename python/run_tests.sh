#!/bin/bash
# Run Python-native tests for LINE solver
# Usage: ./run_native.sh [options]
#   No arguments:      Run all native tests
#   -h, --help:        Show this help message
#   -q, --quick:       Run only quick tests (exclude slow native tests)
#   -f, --file FILE:   Run specific test file
#   -k, --keyword KW:  Run tests matching keyword
#   -x, --stop:        Stop on first failure
#   -v, --verbose:     Verbose output

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
PYTEST_ARGS="-v"
TEST_FILTER=""
SHOW_HELP=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -h|--help)
            SHOW_HELP=1
            shift
            ;;
        -q|--quick)
            TEST_FILTER="not slow"
            shift
            ;;
        -f|--file)
            TEST_FILE="$2"
            shift 2
            ;;
        -k|--keyword)
            TEST_FILTER="$2"
            shift 2
            ;;
        -x|--stop)
            PYTEST_ARGS="$PYTEST_ARGS -x"
            shift
            ;;
        -v|--verbose)
            PYTEST_ARGS="$PYTEST_ARGS -vv"
            shift
            ;;
        *)
            echo "Unknown option: $1"
            SHOW_HELP=1
            shift
            ;;
    esac
done

if [ $SHOW_HELP -eq 1 ]; then
    echo "Usage: $0 [options]"
    echo ""
    echo "Run Python-native tests (same tests as wrapper but with native Python solvers)"
    echo ""
    echo "Options:"
    echo "  -h, --help      Show this help message"
    echo "  -q, --quick     Run only quick tests (exclude slow tests)"
    echo "  -f, --file FILE Run specific test file (relative to tests/)"
    echo "  -k, --keyword K Run tests matching keyword expression"
    echo "  -x, --stop      Stop on first failure"
    echo "  -v, --verbose   Verbose output (extra verbosity)"
    echo ""
    echo "Examples:"
    echo "  $0                              # Run all native tests"
    echo "  $0 --quick                      # Run quick tests only"
    echo "  $0 -f test_all_examples.py      # Run specific file"
    echo "  $0 -k 'gettingstarted'          # Run tests matching 'gettingstarted'"
    exit 0
fi

echo "=========================================="
echo "LINE Solver Python-Native Test Runner"
echo "=========================================="
echo ""

# Force native Python solver implementations
export LINE_SOLVER_LANG=python
# Enable strict mode to disable fallback from native to wrapper (see _backend.py)
export LINE_NATIVE=true
export LINE_NATIVE_STRICT=true
# Fix seed for simulation reproducibility
export LINE_SEED=23000

# Check if Poetry is available
if ! command -v poetry &> /dev/null; then
    echo -e "${RED}Error: poetry not found${NC}"
    exit 1
fi

# Install/update dependencies
echo -e "${BLUE}Installing dependencies...${NC}"
poetry install --no-interaction 2>&1 | grep -E "(Installing|Skipping|up to date)" || true
echo ""

# Determine what to run
# Native tests now run the same test files as wrapper tests, but with LINE_SOLVER_LANG=python
if [ -n "$TEST_FILE" ]; then
    echo -e "${YELLOW}Running native test file: $TEST_FILE${NC}"
    TEST_PATH="tests/$TEST_FILE"
elif [ -n "$TEST_FILTER" ]; then
    echo -e "${YELLOW}Running native tests matching: $TEST_FILTER${NC}"
    TEST_PATH="tests/"
else
    echo -e "${YELLOW}Running all native tests...${NC}"
    TEST_PATH="tests/"
fi

echo ""

# Run tests - disable set -e to capture exit code properly
set +e
if [ -n "$TEST_FILTER" ]; then
    poetry run python3 -m pytest $TEST_PATH -k "$TEST_FILTER" $PYTEST_ARGS
else
    poetry run python3 -m pytest $TEST_PATH $PYTEST_ARGS
fi
TEST_EXIT_CODE=$?
set -e

echo ""
echo "=========================================="
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}All native tests passed!${NC}"
else
    echo -e "${RED}Some native tests failed (exit code: $TEST_EXIT_CODE)${NC}"
fi
echo "=========================================="

# Clean up environment variables
unset LINE_SOLVER_LANG
unset LINE_NATIVE
unset LINE_NATIVE_STRICT
unset LINE_SEED

exit $TEST_EXIT_CODE
