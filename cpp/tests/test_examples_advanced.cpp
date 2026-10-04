/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 *
 * Regression over the SHIPPED ADVANCED examples (cpp/examples/advanced/), the
 * C++ twin of jline.examples.ExampleRunnerAdvancedTest. See
 * test_examples_basic.cpp for the split, and example_suite.h for what is
 * asserted.
 */
#include "example_suite.h"

/** 83 registrations across 16 source files as of 2026-09-12. See the basic twin. */
TEST_CASE("examples/advanced: every shipped example runs and produces finite numbers") {
    line_tests::run_examples_under("advanced/", 83);
}
