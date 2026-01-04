#!/bin/bash
#
# LINE REST API Demo Script
# Demonstrates the capabilities of the LINE REST API
#

set -e

# Default port
PORT=9090
BASE_URL=""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --port|-p)
            PORT="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: ./demo.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -p, --port PORT    Server port (default: 9090)"
            echo "  -h, --help         Show this help"
            echo ""
            echo "Make sure the server is running first with: ./run.sh --port PORT"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

BASE_URL="http://localhost:${PORT}/api/v1"

# Helper function to print section headers
section() {
    echo ""
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${BLUE}  $1${NC}"
    echo -e "${BLUE}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
}

# Helper function for API calls
api_call() {
    local method=$1
    local endpoint=$2
    local data=$3

    echo -e "${YELLOW}>>> $method ${BASE_URL}${endpoint}${NC}"
    if [ -n "$data" ]; then
        echo -e "${YELLOW}Request body:${NC}"
        echo "$data" | jq . 2>/dev/null || echo "$data"
    fi
    echo ""

    if [ "$method" = "GET" ]; then
        response=$(curl -s -X GET "${BASE_URL}${endpoint}" -H "Content-Type: application/json")
    elif [ "$method" = "POST" ]; then
        response=$(curl -s -X POST "${BASE_URL}${endpoint}" -H "Content-Type: application/json" -d "$data")
    elif [ "$method" = "DELETE" ]; then
        response=$(curl -s -X DELETE "${BASE_URL}${endpoint}" -H "Content-Type: application/json")
    fi

    echo -e "${GREEN}Response:${NC}"
    echo "$response" | jq . 2>/dev/null || echo "$response"
    echo ""
}

# Check if server is running
check_server() {
    echo -e "${YELLOW}Checking if server is running on port ${PORT}...${NC}"
    if curl -s --connect-timeout 2 "${BASE_URL}/health" > /dev/null 2>&1; then
        echo -e "${GREEN}Server is running!${NC}"
    else
        echo -e "${RED}Error: Server is not running on port ${PORT}${NC}"
        echo "Start the server first with: ./run.sh --port ${PORT}"
        exit 1
    fi
}

# Create a sample JSIMG model (Closed Queueing Network - ideal for MVA)
# This is a closed network with 5 customers circulating between a Delay station and a Queue
create_sample_model() {
    cat << 'JSIMG_MODEL'
<?xml version="1.0" encoding="ISO-8859-1" standalone="no"?>
<archive xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" name="closed_1class_2stat.jsimg"
    timestamp="Tue Sep 19 03:05:12 BST 2017" xsi:noNamespaceSchemaLocation="Archive.xsd">
    <sim disableStatisticStop="false" logDecimalSeparator="." logDelimiter="," logPath=""
        logReplaceMode="0" maxSamples="1000000" name="closed_1class_2stat.jsimg" polling="1.0"
        xsi:noNamespaceSchemaLocation="SIMmodeldefinition.xsd">
        <userClass customers="5" name="Class1" priority="0" referenceSource="Delay 1" type="closed" />
        <node name="Delay 1">
            <section className="Queue">
                <parameter classPath="java.lang.Integer" name="size">
                    <value>-1</value>
                </parameter>
                <parameter array="true" classPath="java.lang.String" name="dropStrategies">
                    <refClass>Class1</refClass>
                    <subParameter classPath="java.lang.String" name="dropStrategy">
                        <value>drop</value>
                    </subParameter>
                </parameter>
                <parameter classPath="jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy"
                    name="FCFSstrategy" />
                <parameter array="true" classPath="jmt.engine.NetStrategies.QueuePutStrategy"
                    name="QueuePutStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy"
                        name="TailStrategy" />
                </parameter>
            </section>
            <section className="Delay">
                <parameter array="true" classPath="jmt.engine.NetStrategies.ServiceStrategy"
                    name="ServiceStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy"
                        name="ServiceTimeStrategy">
                        <subParameter classPath="jmt.engine.random.Exponential" name="Exponential" />
                        <subParameter classPath="jmt.engine.random.ExponentialPar" name="distrPar">
                            <subParameter classPath="java.lang.Double" name="lambda">
                                <value>1.0</value>
                            </subParameter>
                        </subParameter>
                    </subParameter>
                </parameter>
            </section>
            <section className="Router">
                <parameter array="true" classPath="jmt.engine.NetStrategies.RoutingStrategy"
                    name="RoutingStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy"
                        name="Random" />
                </parameter>
            </section>
        </node>
        <node name="Queue 1">
            <section className="Queue">
                <parameter classPath="java.lang.Integer" name="size">
                    <value>-1</value>
                </parameter>
                <parameter array="true" classPath="java.lang.String" name="dropStrategies">
                    <refClass>Class1</refClass>
                    <subParameter classPath="java.lang.String" name="dropStrategy">
                        <value>drop</value>
                    </subParameter>
                </parameter>
                <parameter classPath="jmt.engine.NetStrategies.QueueGetStrategies.FCFSstrategy"
                    name="FCFSstrategy" />
                <parameter array="true" classPath="jmt.engine.NetStrategies.QueuePutStrategy"
                    name="QueuePutStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.QueuePutStrategies.TailStrategy"
                        name="TailStrategy" />
                </parameter>
            </section>
            <section className="PSServer">
                <parameter classPath="java.lang.Integer" name="maxJobs">
                    <value>1</value>
                </parameter>
                <parameter array="true" classPath="java.lang.Integer" name="numberOfVisits">
                    <refClass>Class1</refClass>
                    <subParameter classPath="java.lang.Integer" name="numberOfVisits">
                        <value>1</value>
                    </subParameter>
                </parameter>
                <parameter array="true" classPath="jmt.engine.NetStrategies.ServiceStrategy"
                    name="ServiceStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.ServiceStrategies.ServiceTimeStrategy"
                        name="ServiceTimeStrategy">
                        <subParameter classPath="jmt.engine.random.Exponential" name="Exponential" />
                        <subParameter classPath="jmt.engine.random.ExponentialPar" name="distrPar">
                            <subParameter classPath="java.lang.Double" name="lambda">
                                <value>1.0</value>
                            </subParameter>
                        </subParameter>
                    </subParameter>
                </parameter>
                <parameter array="true" classPath="jmt.engine.NetStrategies.PSStrategy"
                    name="PSStrategy">
                    <refClass>Class1</refClass>
                    <subParameter classPath="jmt.engine.NetStrategies.PSStrategies.EPSStrategy"
                        name="EPSStrategy" />
                </parameter>
                <parameter array="true" classPath="java.lang.Double" name="serviceWeights">
                    <refClass>Class1</refClass>
                    <subParameter classPath="java.lang.Double" name="serviceWeight">
                        <value>1.0</value>
                    </subParameter>
                </parameter>
            </section>
            <section className="Router">
                <parameter array="true" classPath="jmt.engine.NetStrategies.RoutingStrategy"
                    name="RoutingStrategy">
                    <refClass>Class1</refClass>
                    <subParameter
                        classPath="jmt.engine.NetStrategies.RoutingStrategies.RandomStrategy"
                        name="Random" />
                </parameter>
            </section>
        </node>
        <connection source="Delay 1" target="Queue 1" />
        <connection source="Queue 1" target="Delay 1" />
        <preload>
            <stationPopulations stationName="Delay 1">
                <classPopulation population="5" refClass="Class1" />
            </stationPopulations>
        </preload>
    </sim>
</archive>
JSIMG_MODEL
}

# Escape JSON content
escape_json() {
    echo "$1" | python3 -c 'import json,sys; print(json.dumps(sys.stdin.read()))'
}

#------------------------------------------------------------------------------
# MAIN DEMO
#------------------------------------------------------------------------------

echo ""
echo -e "${GREEN}╔═══════════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${GREEN}║                     LINE REST API DEMO                                    ║${NC}"
echo -e "${GREEN}║                                                                           ║${NC}"
echo -e "${GREEN}║  This script demonstrates the capabilities of the LINE REST API          ║${NC}"
echo -e "${GREEN}║  for solving queueing network models.                                    ║${NC}"
echo -e "${GREEN}╚═══════════════════════════════════════════════════════════════════════════╝${NC}"
echo ""
echo -e "${YELLOW}Using port: ${PORT}${NC}"
echo ""

check_server

# Helper function to pause between sections
pause() {
    echo ""
    echo -e "${YELLOW}Press Enter to continue...${NC}"
    read -r
}

#------------------------------------------------------------------------------
section "1. HEALTH & STATUS ENDPOINTS"
#------------------------------------------------------------------------------

echo "These endpoints are useful for monitoring and Kubernetes health checks."
echo ""

echo -e "${YELLOW}1.1 Health Check (Liveness Probe)${NC}"
api_call "GET" "/health"

echo -e "${YELLOW}1.2 Readiness Check${NC}"
api_call "GET" "/ready"

echo -e "${YELLOW}1.3 Server Info${NC}"
api_call "GET" "/info"

pause

#------------------------------------------------------------------------------
section "2. SOLVER INFORMATION"
#------------------------------------------------------------------------------

echo "LINE supports multiple solvers for different types of analysis."
echo ""

echo -e "${YELLOW}2.1 List All Available Solvers${NC}"
api_call "GET" "/solvers"

echo -e "${YELLOW}2.2 Get Details for MVA Solver${NC}"
api_call "GET" "/solvers/mva"

pause

#------------------------------------------------------------------------------
section "3. MODEL VALIDATION"
#------------------------------------------------------------------------------

echo "Validate models before solving to check for syntax and semantic errors."
echo ""

# Get the sample model
MODEL_CONTENT=$(create_sample_model)
ESCAPED_MODEL=$(escape_json "$MODEL_CONTENT")

echo -e "${YELLOW}3.1 Validate a JSIMG Model${NC}"
VALIDATE_REQUEST=$(cat <<EOF
{
    "format": "jsimg",
    "content": $ESCAPED_MODEL
}
EOF
)
api_call "POST" "/models/validate" "$VALIDATE_REQUEST"

pause

#------------------------------------------------------------------------------
section "4. SOLVING MODELS"
#------------------------------------------------------------------------------

echo "The core functionality: solve queueing network models using various algorithms."
echo ""

echo -e "${YELLOW}4.1 Solve with MVA (Mean Value Analysis)${NC}"
echo "MVA is an analytical solver that computes average performance metrics."
echo ""

SOLVE_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "mva",
    "analysis": "avg"
}
EOF
)
api_call "POST" "/models/solve" "$SOLVE_REQUEST"

echo -e "${YELLOW}4.2 Solve with CTMC (Continuous-Time Markov Chain)${NC}"
echo "CTMC provides exact state-space analysis."
echo ""

SOLVE_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "ctmc",
    "analysis": "avg"
}
EOF
)
api_call "POST" "/models/solve" "$SOLVE_REQUEST"

echo -e "${YELLOW}4.3 Solve with Fluid Solver${NC}"
echo "The fluid solver uses ODEs for mean-field approximation."
echo ""

SOLVE_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "fluid",
    "analysis": "avg"
}
EOF
)
api_call "POST" "/models/solve" "$SOLVE_REQUEST"

pause

#------------------------------------------------------------------------------
section "5. ASYNCHRONOUS JOB SUBMISSION"
#------------------------------------------------------------------------------

echo "For long-running simulations, use async job submission."
echo ""

echo -e "${YELLOW}5.1 Submit Async Job${NC}"
ASYNC_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "mva",
    "analysis": "avg"
}
EOF
)
response=$(curl -s -X POST "${BASE_URL}/models/solve/async" \
    -H "Content-Type: application/json" \
    -d "$ASYNC_REQUEST")
echo -e "${GREEN}Response:${NC}"
echo "$response" | jq .
JOB_ID=$(echo "$response" | jq -r '.jobId')
echo ""

if [ "$JOB_ID" != "null" ] && [ -n "$JOB_ID" ]; then
    echo -e "${YELLOW}5.2 Check Job Status${NC}"
    sleep 1
    api_call "GET" "/jobs/${JOB_ID}"

    echo -e "${YELLOW}5.3 List All Jobs${NC}"
    api_call "GET" "/jobs"
fi

pause

#------------------------------------------------------------------------------
section "6. WHAT-IF ANALYSIS"
#------------------------------------------------------------------------------

echo "Perform parameter sweeps to analyze system behavior under varying conditions."
echo ""

echo -e "${YELLOW}6.1 What-If Analysis: Vary Service Rate${NC}"
WHATIF_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "mva",
    "analysis": "avg",
    "parameter": {
        "type": "service_rate",
        "station": "Queue 1",
        "class": "Class1"
    },
    "values": [0.5, 1.0, 1.5, 2.0]
}
EOF
)
api_call "POST" "/analysis/whatif" "$WHATIF_REQUEST"

pause

#------------------------------------------------------------------------------
section "7. BOTTLENECK ANALYSIS"
#------------------------------------------------------------------------------

echo "Identify performance bottlenecks in your queueing network."
echo ""

echo -e "${YELLOW}7.1 Detect Bottlenecks${NC}"
BOTTLENECK_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "mva"
}
EOF
)
api_call "POST" "/analysis/bottleneck" "$BOTTLENECK_REQUEST"

pause

#------------------------------------------------------------------------------
section "8. SENSITIVITY ANALYSIS"
#------------------------------------------------------------------------------

echo "Analyze how sensitive performance metrics are to parameter changes."
echo ""

echo -e "${YELLOW}8.1 Sensitivity Analysis${NC}"
SENSITIVITY_REQUEST=$(cat <<EOF
{
    "model": {
        "format": "jsimg",
        "content": $ESCAPED_MODEL
    },
    "solver": "mva",
    "parameters": [
        {
            "type": "service_rate",
            "station": "Queue 1",
            "class": "Class1"
        }
    ],
    "delta": 0.1
}
EOF
)
api_call "POST" "/analysis/sensitivity" "$SENSITIVITY_REQUEST"

pause

#------------------------------------------------------------------------------
section "DEMO COMPLETE"
#------------------------------------------------------------------------------

echo -e "${GREEN}The LINE REST API provides a comprehensive interface for:${NC}"
echo ""
echo "  - Health monitoring and Kubernetes integration"
echo "  - Multiple solver algorithms (MVA, CTMC, Fluid, JMT, etc.)"
echo "  - Synchronous and asynchronous model solving"
echo "  - Model validation and format conversion"
echo "  - What-if parameter sweeps"
echo "  - Bottleneck detection"
echo "  - Sensitivity analysis"
echo ""
echo -e "${BLUE}For more information, see the API documentation or source code.${NC}"
echo ""
