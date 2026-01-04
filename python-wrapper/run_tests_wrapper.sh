#!/bin/bash
# Run Python wrapper (JAR) tests for LINE solver
# Usage: ./run_wrapper.sh [options]
#   No arguments:      Run all wrapper tests
#   -h, --help:        Show this help message
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
    echo "Run Python wrapper (JAR) tests (tests/*.py, excluding tests/native/)"
    echo ""
    echo "Options:"
    echo "  -h, --help      Show this help message"
    echo "  -f, --file FILE Run specific test file (relative to tests/)"
    echo "  -k, --keyword K Run tests matching keyword expression"
    echo "  -x, --stop      Stop on first failure"
    echo "  -v, --verbose   Verbose output (extra verbosity)"
    echo ""
    echo "Examples:"
    echo "  $0                              # Run all wrapper tests"
    echo "  $0 -f test_all_examples.py      # Run specific file"
    echo "  $0 -k 'gettingstarted'          # Run tests matching 'gettingstarted'"
    exit 0
fi

echo "=========================================="
echo "LINE Solver Python-Wrapper Test Runner"
echo "=========================================="
echo ""

# Force Java/JAR solver implementations
export LINE_SOLVER_LANG=java

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
if [ -n "$TEST_FILE" ]; then
    echo -e "${YELLOW}Running wrapper test file: $TEST_FILE${NC}"
    TEST_PATH="tests/$TEST_FILE"
elif [ -n "$TEST_FILTER" ]; then
    echo -e "${YELLOW}Running wrapper tests matching: $TEST_FILTER${NC}"
    TEST_PATH="tests/"
    PYTEST_ARGS="$PYTEST_ARGS --ignore=tests/native/"
else
    echo -e "${YELLOW}Running all wrapper tests...${NC}"
    TEST_PATH="tests/"
    PYTEST_ARGS="$PYTEST_ARGS --ignore=tests/native/"
fi

echo ""

# Run tests - use eval to properly handle quoted arguments
if [ -n "$TEST_FILTER" ]; then
    poetry run python3 -m pytest $TEST_PATH -k "$TEST_FILTER" $PYTEST_ARGS
else
    poetry run python3 -m pytest $TEST_PATH $PYTEST_ARGS
fi

TEST_EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}All wrapper tests passed!${NC}"
else
    echo -e "${RED}Some wrapper tests failed (exit code: $TEST_EXIT_CODE)${NC}"
fi
echo "=========================================="

# Clean up environment variable
unset LINE_SOLVER_LANG

exit $TEST_EXIT_CODE
