#!/bin/bash
# Runs all LINE Java benchmarks
# Output format matches MATLAB allBench

mvn compile exec:java -Pbench -DtmpBuild=true 2>/dev/null
