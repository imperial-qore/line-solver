#!/bin/bash
# Run MATLAB tests for LINE solver

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Create logs directory if it doesn't exist
LOG_DIR="$SCRIPT_DIR/logs"
mkdir -p "$LOG_DIR"

# Generate timestamped log filename
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="$LOG_DIR/test_run_${TIMESTAMP}.txt"

echo "=========================================="
echo "LINE Solver MATLAB Test Runner"
echo "=========================================="
echo "Log file: $LOG_FILE"
echo ""

# Run MATLAB tests with output displayed and saved to log file
matlab -batch "lineStart; warning off; addpath(genpath('../../line-test.git')); allTests" 2>&1 | tee "$LOG_FILE"
TEST_EXIT_CODE=${PIPESTATUS[0]}

echo ""
echo "=========================================="
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo "All tests completed!"
else
    echo "Tests failed (exit code: $TEST_EXIT_CODE)"
fi
echo "Log saved to: $LOG_FILE"
echo "=========================================="

exit $TEST_EXIT_CODE
