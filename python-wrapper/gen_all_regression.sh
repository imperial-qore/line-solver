#!/bin/bash
# Generate regression baselines for Python tests using JAR wrapper
# Usage: ./gen_all_regression.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo "=========================================="
echo "LINE Solver Regression Baseline Generator"
echo "(JAR Wrapper - Requires Java)"
echo "=========================================="
echo ""

# Check if Poetry is available
if ! command -v poetry &> /dev/null; then
    echo -e "${RED}Error: poetry not found${NC}"
    exit 1
fi

# Install/update dependencies
echo -e "${BLUE}Installing dependencies...${NC}"
poetry install --no-interaction 2>&1 | grep -E "(Installing|Skipping|up to date)" || true
echo ""

echo -e "${YELLOW}Generating regression baselines using JAR wrapper...${NC}"
echo ""

poetry run python3 tests/gen_all_regression.py

GEN_EXIT_CODE=$?

echo ""
echo "=========================================="
if [ $GEN_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}Regression baselines generated successfully!${NC}"
else
    echo -e "${RED}Regression generation failed (exit code: $GEN_EXIT_CODE)${NC}"
fi
echo "=========================================="

exit $GEN_EXIT_CODE
