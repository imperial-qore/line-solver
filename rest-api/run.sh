#!/bin/bash
#
# LINE REST API Server - Build and Run Script
#

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Check for jline.jar
if [ ! -f "common/jline.jar" ]; then
    echo "Error: common/jline.jar not found"
    echo "Please copy jline.jar to the common/ directory"
    exit 1
fi

# Parse arguments
BUILD=true
PORT=9090
API_KEYS=""
RATE_LIMIT=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --skip-build)
            BUILD=false
            shift
            ;;
        --port)
            PORT="$2"
            shift 2
            ;;
        --api-keys)
            API_KEYS="--api-keys $2"
            shift 2
            ;;
        --rate-limit)
            RATE_LIMIT="--rate-limit $2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: ./run.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skip-build       Skip Maven build"
            echo "  --port PORT        Server port (default: 9090)"
            echo "  --api-keys KEYS    Comma-separated API keys"
            echo "  --rate-limit N     Max requests per minute"
            echo "  -h, --help         Show this help"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Build
if [ "$BUILD" = true ]; then
    echo "Building..."
    mvn clean package -DskipTests -q
fi

# Check JAR exists
if [ ! -f "target/line-rest.jar" ]; then
    echo "Error: target/line-rest.jar not found. Run without --skip-build"
    exit 1
fi

# Run
echo "Starting LINE REST API on port $PORT..."
java -cp "target/line-rest.jar:common/jline.jar" \
    jline.rest.LineRestServer \
    --port "$PORT" $API_KEYS $RATE_LIMIT
